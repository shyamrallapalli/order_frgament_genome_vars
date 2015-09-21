#encoding: utf-8

require 'csv'
require 'bio'
require 'bio-samtools'
require 'yaml'
require 'fileutils'

# input 1 = fasta file
# input 2 = vcf file

# test inputs
# input T1 = ordered fasta file
# input T2 = ordered vcf file
### Read sequence fasta file and store sequences in a hash

sequences = Hash.new {|h,k| h[k] = {} }
file = Bio::FastaFormat.open(ARGV[0])
file.each do |seq|
	sequences[seq.entry_id][:seq] = seq.entry
	sequences[seq.entry_id][:len] = seq.length
end

in_vars = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash of variants from background vcf file
File.open(ARGV[1], 'r').each do |line|
   next if line =~ /^#/
   v = Bio::DB::Vcf.new(line)
   chrom = v.chrom
   pos = v.pos
   info = v.info
   if info["HET"] == "1" or info["AF"] == "0.5"
      in_vars[chrom][:het][pos] = 1
   elsif info["HOM"] == "1" or info["AF"] == "1.0"
      in_vars[chrom][:hom][pos] = 1
   end
end

adj = 1 # ratio adjustments
ratios = []
sequences.keys.each { | id |
	if in_vars[id.to_s].has_key?(:het)
		sequences[id][:het] = in_vars[id.to_s][:het].keys.length + adj
	else
		sequences[id][:het] = adj
	end
	if in_vars[id.to_s].has_key?(:hom)
		sequences[id][:hom] = in_vars[id.to_s][:hom].keys.length + adj
	else
		sequences[id][:hom] = adj
	end
	ratio = (sequences[id][:hom].to_f/sequences[id][:het].to_f)/(sequences[id][:len].to_f/1000000.0) ## normalized ratio per megabase (MB)
	sequences[id][:ratio] = ratio
	ratios << ratio
}

maximum = ratios.max
cutoff = 0.5 # 0.5, 0.2, 0.1, 0.05 times the maximum
minimum = maximum * cutoff

selected = []
sequences.keys.each { | id |
	if sequences[id][:ratio] >= minimum
		selected << id
	end
}


