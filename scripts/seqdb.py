from bidict import *
import os.path, os
import re, sys, util

class SeqDB():
	
	def __init__(self, fn):

		# Initialize attributes
		self.fn = fn # db filename
		self.db = bidict({}) # otu <-> seq
		self.size = {} # otu -> size
		
		# Load sequence database
		#self = self.load_db()
	
	
	def load_db(self):
		# Load existing SeqDB (if exists)
		if os.path.exists(self.fn):
			for tag, seq in util.iter_fst(self.fn):
				otu, size = re.search('>(.*);size=(\d+)', tag).groups()
				self.db[int(otu)] = seq
				self.size[int(otu)] = int(size)
		return self
	
	def add_seq(self, seq, size=1):
		# Add new sequence to SeqDB
		if seq not in ~self.db:
			# If SeqDB is empty, set OTU = 1
			if len(self.db) == 0:
				otu = 1
			# Otherwise, increment to get next OTU
			else:
				otu = max(self.db) + 1
			self.db[otu] = seq
			self.size[otu] = size
		else:
			otu = self.db[:seq]
			self.size[otu] += size
		return otu
	
	
	def merge_db(self, x, keep=0):
		# Merge SeqDB with another database
		# If keep == 0, use ids in self
		if keep == 0:
			for otu in x.db:
				seq = x.db[otu]
				size = x.size[otu]
				self.add_seq(seq, size=size)
			return self
		# If keep == 1, use ids in x
		elif keep == 1:
			for otu in self.db:
				seq = self.db[otu]
				size = self.size[otu]
				x.add_seq(seq, size=size)
			return x
	
	
	def get_seq(self, otu):
		# Get sequence associated with OTU id
		try:
			seq = self.db[otu]
			return seq
		except:
			quit('Error: otu "%s "not in database' %(otu))
	
	
	def get_otu(self, seq, size=1):
		# Get OTU id associated with sequence
		# If sequence in SeqDB, get OTU id
		if seq in ~self.db:
			otu = self.db[:seq]
		# Otherwise, create new SeqDB entry
		else:
			otu = self.add_seq(seq, size=size)
		# Return OTU id
		return otu
	
		
	def trim_db(self, l, keep_all=False):
		
		# ERROR : Identical sequences will double size
		
		# Trim sequences in SeqDB to length l
		for otu in self.db:
			seq = self.db[otu]
			size = self.size[otu]
			new_seq = seq[:l]
			self.add_seq(new_seq, size=size)
			# Remove other sequences from SeqDB (unless keep_all==True)
			if keep_all == False:
				if len(seq) != l:
					del self.db[otu]
		return self
	
	
	def validate(self, fn):
		# Load file as SeqDB
		x = SeqDB(fn)
                # Load sequence database
		x = x.load_db()
		if self.db == x.db:
			return True
		else:
			return False
	
	
	def write(self, out_fn=None, fmt='fasta',  overwrite=False, sort=True):
		
		# If no outfile, use infile
		if out_fn is None:
			out_fn = self.fn
			while overwrite == False and os.path.exists(out_fn):
				out_fn += '.bk'
		
		# Write SeqDB to tempfile
		tmp_fn = '%s.tmp' %(out_fn)
		tmp = open(tmp_fn, 'w')
		for otu in sorted(self.db, key=lambda x: self.size[x], reverse=True):
			seq = self.db[otu]
			size = self.size[otu]
			if fmt == 'fasta':
				tmp.write('>%d;size=%d\n%s\n' %(otu, size, seq))
			elif fmt == 'tab':
				tmp.write('%d\t%s\n' %(otu, seq))
			else:
				quit('Unrecognized format option')
		tmp.close()
		
		# Validate and move to final destination
		if self.validate(tmp_fn):
			cmd = 'mv %s %s' %(tmp_fn, out_fn)
			os.system(cmd)
		return self
	
	
	def map_db(self, x, reverse=False):
		# Map OTUs in self to OTUs in another SeqDB
		if reverse == False:
			m = {}
			for seq in ~self.db:
				otu1 = self.get_otu(seq)
				otu2 = x.get_otu(seq)
				m[otu1] = otu2
		# Or map OTUs in another SeqDB to self
		elif reverse == True:
			m = {}
			for seq in ~self.x:
				otu1 = x.get_otu(seq)
				otu2 = self.get_otu(seq)
				m[otu1] = otu2
		# Otherwise, throw error
		else:
			quit('error: invalid argument in map_to_db()')
		return m
	
