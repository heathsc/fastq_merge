use std::{
    fmt, fs,
    io::{self, BufRead, BufReader, Error, ErrorKind, Read},
    path::{Path, PathBuf},
};

use compress_io::compress::CompressIo;

/// Shared configuration data set from the command line
pub struct Config {
    report_file: Option<String>,
    inputs: Vec<InputFile>,
    summary_files: Vec<SummaryFile>,
    summary_file_name: String,
}

impl Config {
    /// Creates a 'Config' from the components.
    /// We will build up a vector of fastq files in inputs and a vector of sequence summary files in summary_files.
    ///
    ///
    /// The user supplies either an optional list of input names and an optional list of summary file names.  The
    /// input names can either be individual fastq files or can be directories containing fastq files. Optionally
    /// the names of the report and summary output files can also be supplied.
    ///
    /// If an input name is a file then this is added to inputs.  If the input name is a directory then
    /// it is checked to see if it looks like a guppy run directory in containing a file called "sequencing_summary.txt"
    /// and a subdirectory named "pass".  In this case then any fastq files in pass are added to inputs and
    /// the sequencing_summary.txt file is added to summary_files, otherwise an fastq files in the original directory
    /// are added to inputs. If no fastq files are found for an input directory this is an error and processing will stop.
    ///
    /// We create a 'Vec' of 'InputFile' from the supplied input names as described above. If no names are supplied
    /// inputs will be a 'Vec' with a single 'None' entry, otherwise inputs will be a 'Vec' with one or more 'Some<Box<str>>'
    /// with the names of the potential fastq files. When processing we can simply iterate over inputs without
    /// worrying about whether this is a regular file, a directory or stdin, as these details have already been handled.
    pub fn new(
        inputs: Vec<String>,
        summary_files: Vec<String>,
        summary_file_name: String,
        report_file: Option<String>,
    ) -> io::Result<Self> {
        // Add and supplied summary file names to summary_files vector.  The files are checked at this stage to make sure that
        // they exist, and that we can read the first line (the header) from them.  The contents of the header line are not
        // checked at this stage.	The summary_vector vector can be added to when the inputs are considered if a guppy run directory
        // is found.
        let mut v = Vec::with_capacity(summary_files.len());
        for s in summary_files.into_iter() {
            match SummaryFile::new(s)? {
                Some(x) => v.push(x),
                None => warn!("Summary file is empty, will be ignored"),
            }
        }
        let mut summary_files = v;

        // Check input names.  If no name supplied create a vector with a single None entry otherwise iterate over each supplied input name
        let inputs = if inputs.is_empty() {
            vec![InputFile::none()]
        } else {
            let mut in_files = Vec::new();
            for f in inputs.into_iter() {
                add_input(f, &mut in_files, &mut summary_files)?;
            }
            in_files
        };

        // Create 'Config'
        Ok(Self {
            inputs,
            summary_files,
            summary_file_name,
            report_file,
        })
    }

    /// Returns a slice with the list of input filenames
    pub fn inputs(&self) -> &[InputFile] {
        &self.inputs
    }

    /// Returns options reference to report file filename
    pub fn report_file(&self) -> Option<&str> {
        self.report_file.as_deref()
    }

    /// Returns a slice with list of summary files
    pub fn summary_files(&self) -> &[SummaryFile] {
        &self.summary_files
    }

    /// Returns output summary file name
    pub fn summary_file_name(&self) -> &str {
        &self.summary_file_name
    }
}

/// Input files from command line
/// If name is 'None' then this represents stdin, otherwise it represents a file
#[derive(Debug)]
pub struct InputFile {
    name: Option<PathBuf>,
}

impl InputFile {
    /// Create new 'InputFile' from a file name as 'Box<str>'
    fn new(name: String) -> Self {
        Self {
            name: Some(PathBuf::from(name)),
        }
    }

    /// Create new 'InputFile' from a file name as PathBuf
    fn from_pathbuf(name: PathBuf) -> Self {
        Self { name: Some(name) }
    }
    /// Create new 'InputFile' representing stdin
    fn none() -> Self {
        Self { name: None }
    }

    /// Get BufReader from 'InputFile'
    /// Returns error if file can not be opened for reading
    pub fn get_bufreader(&self) -> io::Result<BufReader<Box<dyn Read>>> {
        let mut c = &mut CompressIo::new();
        if let Some(p) = self.name.as_deref() {
            c = c.path(p);
        }
        c.bufreader()
    }
}

impl fmt::Display for InputFile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.name {
            Some(s) => write!(f, "{}", s.display()),
            None => write!(f, "<stdin>"),
        }
    }
}

/// Sequence summary input file
/// We store the header line when we first encounter the file.  This is because different files don't always have the same columns
/// so we will need to match the column names when we merge the files to ensure the columns match.
pub struct SummaryFile {
    name: PathBuf,
    header_line: String,
}

impl SummaryFile {
    /// Create new 'SummaryFile' from String.  Returns and error if file cannot be opened
    /// or the first line can not be read.  Returns None if file is empty.
    fn new(name: String) -> io::Result<Option<Self>> {
        Self::from_pathbuf(PathBuf::from(name))
    }

    /// As new() but from PathBuf
    fn from_pathbuf(name: PathBuf) -> io::Result<Option<Self>> {
        let mut rdr = CompressIo::new().path(&name).bufreader()?;
        let mut header_line = String::new();
        Ok(if rdr.read_line(&mut header_line)? == 0 {
            None
        } else {
            Some(Self { name, header_line })
        })
    }

    /// Get BufReader from
    /// Returns error if file can not be opened for reading
    pub fn get_reader(&self) -> io::Result<BufReader<Box<dyn Read>>> {
        CompressIo::new().path(&self.name).bufreader()
    }

    /// Returns stored header line
    pub fn header_line(&self) -> &str {
        &self.header_line
    }
}

impl fmt::Display for SummaryFile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name.display())
    }
}

/// Add input from command line
/// This could be a file (expected to be a fastq file) or a directory.
/// If a directory, this could directly contain a list of fastq files, or it could be
/// the output directory of a guppy basecalling run containing a sequencing_summary.txt file and
/// a directory called "pass" which itself contains fastq files.  Any fastq files will be added to in_files
/// in the order that they are found, and the sequencing_summary file, if found, will be added to summary_files.
/// If no fastq files are found in a supplied directory then this is reported as an error.
fn add_input(
    name: String,
    in_files: &mut Vec<InputFile>,
    summary_files: &mut Vec<SummaryFile>,
) -> io::Result<()> {
    trace!("Checking input: {}", name);
    if fs::metadata(&name)?.is_dir() {
        add_directory(name, in_files, summary_files)?
    } else {
        debug!("Adding input fastq file {}", name);
        in_files.push(InputFile::new(name))
    }
    Ok(())
}

/// Add input files from directory
/// We first check if the file sequencing_summary.txt and sub_directory pass exist in the directory
/// If so we add the sequencing_summary.txt file to summary_files, and add any fastq files found in pass to
/// in_files, otherwise just add any fastq files in the directory to in_files
fn add_directory(
    name: String,
    in_files: &mut Vec<InputFile>,
    summary_files: &mut Vec<SummaryFile>,
) -> io::Result<()> {
    trace!(
        "Checking whether directory looks like a guppy run directory: {}",
        name
    );
    let name = PathBuf::from(name);
    let summ_file = name.join("sequencing_summary.txt");
    trace!("Checking for summary file {}", summ_file.display());
    if summ_file.exists() {
        let pass_dir = name.join("pass");
        trace!(
            "Summary file found.  Checking for pass sub-directory {}",
            pass_dir.display()
        );
        if pass_dir.is_dir() {
            trace!("pass directory found.  Add summary file to list for processing, and check for fastq files in {}", pass_dir.display());
            match SummaryFile::from_pathbuf(summ_file)? {
                Some(s) => {
                    debug!(
                        "Adding discovered summary file {} to list for processing",
                        s
                    );
                    summary_files.push(s);
                }
                None => warn!("Summary file is empty - will be ignored"),
            }
            add_simple_directory(pass_dir, in_files)
        } else {
            trace!(
                "pass directory not present.  Just check for fastq files in {}",
                name.display()
            );
            add_simple_directory(name, in_files)
        }
    } else {
        trace!(
            "Summmary file not present.  Just check for fastq files in {}",
            name.display()
        );
        add_simple_directory(name, in_files)
    }
}

/// Add fastq files from directory
/// Look for files ending in .fastq or .fq with or without a possible compression suffix (i.e, *.fastq.gz).
/// An error is returned either if the directory can not be opened or if no fastq files are found
fn add_simple_directory(name: PathBuf, in_files: &mut Vec<InputFile>) -> io::Result<()> {
    // Track how many fastq files have been found
    let mut n = 0;

    // Iterate over directory entries
    for entry in name.read_dir()? {
        let file = entry?.path();

        // Check whether file looks like a regular or compressed fastq file.  We do this by first looking for a "fastq" or "fq"
        // extension. If not found, remove any file extension and check again.
        if is_fastq(&file) || file.file_stem().map_or(false, is_fastq) {
            debug!("Adding input fastq file {}", file.display());
            in_files.push(InputFile::from_pathbuf(file));
            n += 1;
        }
    }

    // Check that we have found some fastq files
    if n > 0 {
        debug!("Added {} fastq files from directory {}", n, name.display());
        Ok(())
    } else {
        Err(Error::new(
            ErrorKind::Other,
            format!("No fastq files found in directory {}", name.display()),
        ))
    }
}

/// Check if filename ends with a "fastq" or "fq" extension
fn is_fastq<P: AsRef<Path>>(file: P) -> bool {
    let file = file.as_ref();
    let res = file
        .extension()
        .map_or(false, |s| s == "fastq" || s == "fq");
    trace!("is_fastq({}) = {}", file.display(), res);
    res
}
