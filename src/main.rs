use anyhow::{Context, Error};
use chrono::Local;
use clap::Parser;
use indexmap::IndexMap;
use rust_htslib::bam::{
    Format, Header, Read, Reader, Record, Writer, ext::BamRecordExtensions, record::Aux,
    record::Cigar,
};
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Parser, Debug)]
pub struct Args {
    #[arg(short = 'i')]
    input: String,

    #[arg(short = 'd', default_value = "split")]
    output_dir: String,

    #[arg(short = 't', default_value = "BX")]
    tag: String,

    #[arg(short = 'o', default_value = "report.tsv")]
    report_file: String,

    #[arg(short = 's', default_value = "_")]
    separator: String,
}

pub fn validate_args(args: &Args) -> Result<(), Error> {
    // println! {"{:?}", args}
    if args.tag.len() != 2 {
        anyhow::bail!("Tag must be 2 letters!: {}", args.tag.len())
    }
    Ok(())
}

/// Record the number of differences between two read sequences. We start comparison at the latest start coordinate, and walk along both reads until one of them runs out of bases.
///
/// It should be noted that soft- and hard-clipped bases are not compared, and that bases outside the span of the shortest read's mapping coordinates are not compared.
pub fn compare_read_sequences(r1: &Record, r2: &Record) -> usize {
    let mut num_diff = 0;

    let start_pos = std::cmp::max(r1.pos(), r2.pos());

    let mut iter1 = r1.aligned_pairs().map(|s| (s[0], s[1])).peekable();
    let mut iter2 = r2.aligned_pairs().map(|s| (s[0], s[1])).peekable();

    // catch up both iterators to the start position
    while let Some((_read_pos, ref_pos)) = iter1.peek() {
        if *ref_pos < start_pos {
            iter1.next();
        } else {
            break;
        }
    }

    while let Some((_read_pos, ref_pos)) = iter2.peek() {
        if *ref_pos < start_pos {
            iter2.next();
        } else {
            break;
        }
    }

    loop {
        match (iter1.peek(), iter2.peek()) {
            (None, _) => break, // either iterator is exhausted, so break
            (_, None) => break,

            // compare both
            (Some(&(read_a, ref_a)), Some(&(read_b, ref_b))) => {
                // we went out of sync, catch up and register the lag as mismatches
                if ref_a < ref_b {
                    num_diff += 1;
                    iter1.next();
                    continue;
                }

                if ref_b < ref_a {
                    num_diff += 1;
                    iter2.next();
                    continue;
                }

                // we are at the current position
                assert_eq!(
                    ref_a, ref_b,
                    "read positions out of sync in overlap detection"
                );

                if r1.seq()[read_a as usize] != r2.seq()[read_b as usize] {
                    num_diff += 1;
                }

                iter1.next();
                iter2.next();
            }
        }
    }

    num_diff
}

pub fn extract_umi_from_header<'a>(header: &'a str, separator: &str) -> Result<&'a str, Error> {
    let (_rest, past_sep) = header.rsplit_once(separator).with_context(|| {
        format!(
            "failed to get UMI with separator '{}'. Header in question:\n{}",
            separator, header
        )
    })?;

    // check if we have r1/r2 extensions, e.g:
    // <REST_OF_HEADER>_<UMI>/1
    // if we don't remove that, we'll over-stratify read groups.
    if let Some((umi, _mate_info)) = past_sep.rsplit_once("/") {
        Ok(umi)
    } else {
        Ok(past_sep)
    }
}

pub fn get_umi(r: &Record, separator: &str) -> Result<String, Error> {
    unsafe {
        let s = std::str::from_utf8_unchecked(r.qname());
        Ok(String::from(extract_umi_from_header(s, separator)?))
    }
}

pub fn qname_to_string(qname: &[u8]) -> String {
    unsafe { String::from_utf8_unchecked(qname.to_vec()) }
}

pub type TagWriter = (Option<String>, Option<Writer>);

pub fn make_writer_from_file(
    infile: &str,
    tag: &str,
    header: &Header,
    outdir: &str,
    cur_file_name: &mut String,
) -> Result<Writer, Error> {
    let infile = Path::new(infile)
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();

    let (outname, _ext) = infile
        .rsplit_once(".bam")
        .context("Bam file extension not found.")?;
    let writer_fname = format!("{outdir}/{outname}_{tag}.bam");
    let writer = Writer::from_path(&writer_fname, header, Format::Bam)?;
    *cur_file_name = writer_fname.to_string();
    Ok(writer)
}

pub struct VarianceReport {
    diff_count: usize,
    ratio_diff: f32,
    total_count: usize,
    n_diff_pos: usize,
    n_diff_umis: usize,
    n_reads_highest_mapq: usize,
    most_bases_mismatch: usize,
    least_bases_mismatch: usize,
}

pub struct SeqMap {
    inner: IndexMap<u64, SeqEntry>,
    count: usize,
}

impl SeqMap {
    fn intake(&mut self, read: Record) {
        let mut h = std::hash::DefaultHasher::new();
        read.seq().encoded.hash(&mut h);
        let e = self.inner.entry(h.finish()).or_insert(SeqEntry::new());

        e.count += 1;
        self.count += 1;
        e.up_group(read);
    }

    fn new() -> Self {
        let inner = IndexMap::new();
        Self { inner, count: 0 }
    }

    fn contains_indels(&self) -> bool {
        for (_, seq) in self.inner.iter() {
            for r in seq.reads.iter() {
                if read_has_indel(r) {
                    return true;
                }
            }
        }

        false
    }

    // kind of slow right now, need to combine a lot of operations in one loop instead of using
    // several.
    fn tally_deviants(&mut self) -> VarianceReport {
        let mut diff_positions: HashSet<usize> = HashSet::new();
        let mut diff_umis: HashSet<String> = HashSet::new();
        let mut mapq_count: IndexMap<u8, usize> = IndexMap::new();
        let mut least_bases_mismatch = usize::MAX;
        let mut most_bases_mismatch = usize::MIN;

        // sort in reverse by read sequence count, so that 0th kv is highest count
        self.inner
            .sort_by(|_k1, v1, _k2, v2| v2.count.cmp(&v1.count));

        // tally unique UMIs
        self.inner.iter().for_each(|(_k, v)| {
            v.reads.iter().for_each(|r| {
                diff_umis.insert(get_umi(r, "_").expect("Failed to get read UMI"));
                *mapq_count.entry(r.mapq()).or_insert(0) += 1;
            })
        });

        let _maj_count = self.inner.get_index(0).unwrap().1.count;
        let mut diff_count = 0;
        self.inner
            .iter()
            .skip(1)
            .for_each(|(_k, v)| diff_count += v.count);

        // let total_count = (maj_count + diff_count) as usize;
        let total_count = self.count;
        let ratio_diff = diff_count as f32 / (total_count as f32);
        let diff_count = diff_count as usize;

        // let maj_seq = self.inner.get_index(0).unwrap().1.reads[0].seq().as_bytes();
        let maj_seq = &self.inner.get_index(0).unwrap().1.reads[0];

        self.inner.iter().skip(1).for_each(|(_k, v)| {
            let other_seq = &v.reads[0];

            let mismatch_count = compare_read_sequences(maj_seq, other_seq);

            // if we have a read that's longer than the majority sequence, then count the number of
            // extra bases as mismatches.
            // let mut mismatch_count: usize =
            //     (maj_seq.len() as i32 - other_seq.len() as i32).abs() as usize;

            // // note that we only check the bases up to the end of the shortest read.
            // for i in 0..(std::cmp::min(maj_seq.len(), other_seq.len())) {
            //     if maj_seq[i] != other_seq[i] {
            //         diff_positions.insert(i);
            //         mismatch_count += 1;
            //     }
            // }

            least_bases_mismatch = std::cmp::min(mismatch_count, least_bases_mismatch);
            most_bases_mismatch = std::cmp::max(mismatch_count, most_bases_mismatch);
        });

        // no mismatched reads to update mismatch counts, so revert to zero.
        if least_bases_mismatch == usize::MAX {
            least_bases_mismatch = 0;
        }

        let n_diff_pos = diff_positions.len();

        // we subtract one since we always have the majority UMI
        let n_diff_umis = diff_umis.len() - 1;

        // sort to get the highest mapq in 0th idx
        mapq_count.sort_by(|k1, _v1, k2, _v2| k2.cmp(k1));
        let (_mapq, n_reads_highest_mapq) = mapq_count.swap_remove_index(0).unwrap();

        VarianceReport {
            diff_count,
            total_count,
            ratio_diff,
            n_diff_pos,
            n_diff_umis,
            n_reads_highest_mapq,
            most_bases_mismatch,
            least_bases_mismatch,
        }
    }
}

pub struct SeqEntry {
    pub reads: Vec<Record>,
    pub count: i32,
    pub qual_sum: u32,
}

pub fn read_has_indel(r: &Record) -> bool {
    for c in r.cigar().iter() {
        if matches!(c, Cigar::Ins(_)) || matches!(c, Cigar::Del(_)) {
            return true;
        }
    }

    false
}

impl SeqEntry {
    pub fn up_group(&mut self, read: Record) {
        self.reads.push(read);
    }

    pub fn new() -> Self {
        let reads = vec![];
        let count = 0;
        let qual_sum = 0;

        Self {
            reads,
            count,
            qual_sum,
        }
    }
}

pub fn create_file_timestamp(fsuffix: &str) -> String {
    let dt = Local::now();
    let dt = dt.format("%H%M%S_%d-%b-%Y");
    format!("{dt}_{fsuffix}")
}

pub fn write_cluster_report(
    writer: &mut BufWriter<std::fs::File>,
    tag: &str,
    report: VarianceReport,
) -> Result<(), Error> {
    Ok(writeln!(
        writer,
        "{tag}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        report.diff_count,
        report.n_diff_pos,
        report.ratio_diff,
        report.n_diff_umis,
        report.total_count,
        report.n_reads_highest_mapq,
        report.least_bases_mismatch,
        report.most_bases_mismatch,
    )?)
}

/// This assumes an input bam that is sorted by the same tag as used in the runtime args
fn _main() -> Result<(), Error> {
    let args = Args::parse();
    validate_args(&args)?;

    let outdir = Path::new(&args.output_dir);

    if outdir.exists() && outdir.is_file() {
        anyhow::bail!(
            "Output dir {} already exists as a file! Not overwriting. Exiting...",
            args.output_dir
        );
    }

    if !outdir.exists() {
        std::fs::create_dir(outdir)?;
    }

    let report_file = create_file_timestamp(&args.report_file);
    let fhandle = std::fs::File::create_new(report_file)?;
    let mut bufwriter = BufWriter::new(fhandle);
    bufwriter.write_all(b"tag\tnum_diff\tnum_unique_diff_pos\tfrac_diff\tn_diff_umis\ttotal\tn_reads_highest_mapq\tleast_bases_mismatch\tmost_bases_mismatch\n")?;

    let tag = args.tag.as_bytes();
    let mut cur_file_name = "".to_string();

    let mut reader = Reader::from_path(&args.input)?;
    reader.set_threads(num_cpus::get())?;
    let header = Header::from_template(reader.header());

    let mut writer: Writer;
    let mut prev_tag = "".to_string();

    let mut read_store = SeqMap::new();

    let mut n_clusters_used = 0;
    let mut n_clusters_with_indels = 0;

    for r in reader
        .records()
        .map(|rec| rec.expect("Failed to read read!"))
    {
        match r.aux(tag) {
            Ok(value) => {
                if let Aux::String(cur_tag) = value {
                    // new tag doesn't match old one, so we've encountered a new tag group in our
                    // sorted file
                    if cur_tag != prev_tag {
                        if read_store.count >= 1 && prev_tag != "" {
                            // don't use clusters with indels in filtering.
                            // if read_store.contains_indels() {
                            // n_clusters_with_indels += 1;
                            // } else {
                            n_clusters_used += 1;
                            let report = read_store.tally_deviants();
                            write_cluster_report(&mut bufwriter, &prev_tag, report)?;

                            if read_store.count >= 1000 {
                                writer = make_writer_from_file(
                                    &args.input,
                                    &prev_tag,
                                    &header,
                                    &args.output_dir,
                                    &mut cur_file_name,
                                )?;

                                for (_k, v) in read_store.inner.drain(..) {
                                    for r in v.reads {
                                        writer.write(&r)?;
                                    }
                                }
                            }
                        }
                        // }

                        read_store.inner.clear();
                        read_store.count = 0;
                        // println! {"{cur_file_name}"}

                        assert!(read_store.inner.is_empty());

                        prev_tag = cur_tag.to_string();
                        read_store.intake(r);

                        // found another read clustered under the current tag
                    } else {
                        read_store.intake(r);
                    }
                }
            }
            // no tag found for this read.
            Err(_) => continue,
        }
    }

    writer = make_writer_from_file(
        &args.input,
        &prev_tag,
        &header,
        &args.output_dir,
        &mut cur_file_name,
    )?;

    // println! {"last group..."}
    // if read_store.contains_indels() {
    // n_clusters_with_indels += 1;
    // } else {
    n_clusters_used += 1;
    let report = read_store.tally_deviants();
    write_cluster_report(&mut bufwriter, &prev_tag, report)?;
    for (_k, v) in read_store.inner.drain(..) {
        for r in v.reads {
            writer.write(&r)?;
        }
    }
    // }

    println! {"FILE: {}\t CLUSTERS DQed: {n_clusters_with_indels}\t CLUSTERS USED: {n_clusters_used}", args.input};

    Ok(())
}

pub fn main() {
    _main().map_err(|e| eprintln! {"FAIL: {e}"}).unwrap();
}
