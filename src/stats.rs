use std::{
	ops::AddAssign,
	fmt
};

const QUALITY_OFFSET: usize = 33;

/// Collection of basic stats about reads
/// For the moment we are just collecting information
/// on the number of records, the number of bases, 
/// base quality distribution etc. as well as the
/// number of bases with illegal qualities (< QUALITY_OFFSET)
pub struct BasicStats {
	records: usize,
	illegal: usize,
	quality_dist: [usize; 256 - QUALITY_OFFSET],
}

impl Default for BasicStats {
/// Create new empty 'BasicStats'
	fn default() -> Self {
		Self {
			records: 0,
			illegal: 0,
			quality_dist: [0; 256 - QUALITY_OFFSET],
		}
	}	
}

impl AddAssign for BasicStats {
	fn add_assign(&mut self, other: Self) {
		self.records += other.records;
		self.illegal += other.illegal;
		for (q1, q2) in self.quality_dist.iter_mut().zip(other.quality_dist.iter()) {
			*q1 += *q2
		}
	}
}

impl fmt::Display for BasicStats {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { 
		let bases: usize = self.quality_dist.iter().sum();
		let mut sum = 0;
		let med = {
			let mut qx = 0;
			for (q, n) in self.quality_dist.iter().enumerate() {
				sum += n;
				if sum >= bases >> 1 {
					qx = q;
					break
				}
			}
			qx		
		};
		if self.illegal == 0 {
			write!(f, "reads = {}. bases = {}, median_qual = {}", self.records, bases, med)
		} else {
			write!(f, "reads = {}. bases = {}, illegal = {}, median_qual = {}", self.records, bases, self.illegal, med)
		}	
	}
}

impl BasicStats {
/// Updates 'BasicStats' based on quality '&str' from FASTQ record
	pub fn add_record(&mut self, qual: &str) {
		self.records += 1;
		let qual = qual.as_bytes();
		for q in qual.iter().map(|x| *x as usize) {
			if q >= QUALITY_OFFSET {
				self.quality_dist[q - QUALITY_OFFSET] += 1;
			} else {
				self.illegal += 1;
			}
		}
	}
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum FileStatsType { Unique, NonUnique, Duplex }

/// Collection of basic stats about an input file
/// We store a BasicStats for the unique and non-unique reads separately
#[derive(Default)]
pub struct FileStats {
	unique: BasicStats,
	non_unique: BasicStats,
	duplex: BasicStats,
}

impl FileStats {
/// Creates new 'FileStats' with zero counts
	pub fn new() -> Self { Self::default() }
	
/// Updates 'FileStats' based on quality '&str' from FASTQ record.
/// The unique flag indicates whether the record is unique or a duplicate
	pub fn add_record(&mut self, qual: &str, fs_type: FileStatsType) {
		match fs_type {
			FileStatsType::Unique => &mut self.unique,
			FileStatsType::NonUnique => &mut self.non_unique,
			FileStatsType::Duplex => &mut self.duplex,
		}.add_record(qual)
	}
}

impl AddAssign for FileStats {
	fn add_assign(&mut self, other: Self) {
		self.unique += other.unique;
		self.non_unique += other.non_unique;
		self.duplex += other.duplex;
	}
}

impl fmt::Display for FileStats {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { 
		write!(f, "Unique reads: {}, Non-unique reads: {}, Duplex reads: {}", self.unique, self.non_unique, self.duplex)
	}
}

