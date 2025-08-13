use anyhow::{Context, Error};
use chrono::Local;
use clap::Parser;
use indexmap::IndexMap;
use rust_htslib::bam::{Format, Header, Read, Reader, Record, Writer, record::Aux};
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Parser, Debug)]
pub struct Args {
    input: String,
    output_dir: String,
    tag: String,
    report_file: String,
}

pub fn validate_args(args: &Args) -> Result<(), Error> {
    println! {"{:?}", args}
    if args.tag.len() != 2 {
        anyhow::bail!("Tag must be 2 letters!: {}", args.tag.len())
    }
    Ok(())
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

    fn tally_deviants(&mut self) -> VarianceReport {
        let mut diff_positions: HashSet<usize> = HashSet::new();

        // sort in reverse by read sequence count, so that 0th kv is highest count
        self.inner
            .sort_by(|_k1, v1, _k2, v2| v2.count.cmp(&v1.count));

        let _maj_count = self.inner.get_index(0).unwrap().1.count;
        let mut diff_count = 0;
        self.inner
            .iter()
            .skip(1)
            .for_each(|(_k, v)| diff_count += v.count);

        // let total_count = (maj_count + diff_count) as usize;
        let total_count = self.count as usize;
        let ratio_diff = diff_count as f32 / (total_count as f32);
        let diff_count = diff_count as usize;

        let maj_seq = self.inner.get_index(0).unwrap().1.reads[0].seq().as_bytes();

        self.inner.iter().skip(1).for_each(|(_k, v)| {
            let other_seq = v.reads[0].seq().as_bytes();

            // note that we only check the bases up to the end of the shortest read.
            for i in 0..(std::cmp::min(maj_seq.len(), other_seq.len())) {
                if maj_seq[i] != other_seq[i] {
                    diff_positions.insert(i);
                }
            }
        });

        let n_diff_pos = diff_positions.len();

        VarianceReport {
            diff_count,
            total_count,
            ratio_diff,
            n_diff_pos,
        }
    }
}

pub struct SeqEntry {
    pub reads: Vec<Record>,
    pub count: i32,
    pub qual_sum: u32,
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
        "{tag}\t{}\t{}\t{}\t{}",
        report.diff_count, report.n_diff_pos, report.ratio_diff, report.total_count
    )?)
}

/// This assumes an input bam that is sorted by the same tag as used in the runtime args
fn main() -> Result<(), Error> {
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
    bufwriter.write(b"tag\tnum_diff\tnum_unique_diff_pos\tfrac_diff\ttotal\n")?;

    let tag = args.tag.as_bytes();
    let mut cur_file_name = "".to_string();

    let mut reader = Reader::from_path(&args.input)?;
    reader.set_threads(num_cpus::get())?;
    let header = Header::from_template(&reader.header());

    let mut writer: Writer;
    let mut prev_tag = "".to_string();

    let mut read_store = SeqMap::new();

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
                        if read_store.count >= 3 && prev_tag != "" {
                            let report = read_store.tally_deviants();
                            write_cluster_report(&mut bufwriter, &prev_tag, report)?;

                            if read_store.count >= 100 {
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

                        read_store.inner.clear();
                        read_store.count = 0;
                        println! {"{cur_file_name}"}

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

    println! {"last group..."}
    let report = read_store.tally_deviants();
    write_cluster_report(&mut bufwriter, &prev_tag, report)?;
    for (_k, v) in read_store.inner.drain(..) {
        for r in v.reads {
            writer.write(&r)?;
        }
    }

    Ok(())
}
