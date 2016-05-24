import die

# 	Code	Means				Comp	CompCode
# 	{ 'M', "AC",   'K' },		// GT		K
# 	{ 'R', "AG",   'Y' },		// CT		Y
# 	{ 'W', "AT",   'W' },		// AT		W
# 	{ 'S', "CG",   'S' },		// CG		S
# 	{ 'Y', "CT",   'R' },		// AG		R
# 	{ 'K', "GT",   'M' },		// AC		M
# 	{ 'V', "ACG",  'B' },		// CGT		B
# 	{ 'H', "ACT",  'D' },		// AGT		D
# 	{ 'D', "AGT",  'H' },		// ACT		H
# 	{ 'B', "CGT",  'V' },		// ACG		V
# 	{ 'X', "GATC", 'X' },		// ACGT		X
# 	{ 'N', "GATC", 'N' },		// ACGT		N

def MatchLetter(a, b):
	if b == 'A' or b == 'C' or b == 'G' or b == 'T':
		return a == b
	elif b == 'S':
		return a == 'C' or a == 'G'
	elif b == 'Y':
		return a == 'T' or a == 'C'
	elif b == 'K':
		return a == 'G' or a == 'T'
	elif b == 'V':
		return a == 'A' or a == 'C' or a == 'G'
	elif b == 'H':
		return a == 'A' or a == 'C' or a == 'T'
	elif b == 'D':
		return a == 'A' or a == 'G' or a == 'T'
	elif b == 'B':
		return a == 'C' or a == 'G' or a == 'T'
	elif b == 'X' or b == 'N':
		return a == 'G' or a == 'A' or a == 'T' or a == 'C'
	elif b == 'R':
		return a == 'A' or a == 'G'
	elif b == 'M':
		return a == 'A' or a == 'C'
	elif b == 'H':
		return a == 'A' or a == 'C' or a == 'T'
	elif b == 'W':
		return a == 'A' or a == 'T'
	else:
		die.Die("Bad letter in primer '%c'" % b)
	return 0

def MatchPrefix(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	n = PrimerLength
	if L < n:
		n = L
	Diffs = 0
	for i in range(0, n):
		if not MatchLetter(Seq[i], Primer[i]):
			Diffs += 1
	return Diffs

def MatchPos(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	for Pos in range(0, L-PrimerLength):
		if MatchPrefix(Seq[Pos:], Primer) == 0:
			return Pos
	return -1

def Match(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	for Pos in range(0, L-PrimerLength):
		if MatchPrefix(Seq[Pos:], Primer) == 0:
			return Pos
	return -1

def BestMatch(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	BestDiffs = PrimerLength
	BestPos = -1
	for Pos in range(0, L-PrimerLength+1):
		d = MatchPrefix(Seq[Pos:], Primer)
		if d < BestDiffs:
			BestDiffs = d
			BestPos = Pos
	return BestPos, BestDiffs

