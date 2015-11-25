#encoding: utf-8

require 'csv'
require 'bio'
require 'bio-samtools'
require 'yaml'
require 'fileutils'
require "rinruby"

# input 1 = fasta file
# input 2 = vcf file

# test inputs
# input T1 = ordered fasta file
# input T2 = ordered vcf file

# filter parameter are to be read from a file in the current folder
pars = YAML.load_file("./filter_pars.yml")
adj = pars['ratio_adj']  # 1 or 0.1 or 0.01 etc will be added to numertator and denominator
filter = pars['filter']   # cut off filter calculated as percent of maximum ratio

### Read sequence fasta file and store sequences in a hash
sequences = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
file = Bio::FastaFormat.open(ARGV[0])
file.each do |seq|
	sequences[seq.entry_id][:seq] = seq.entry
	sequences[seq.entry_id][:len] = seq.length.to_i
end

### read vcf file and make a hash of variants from vcf file
in_vars = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
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

snpratio = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
### process sequence fragments based on variant density
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
	ratio = 0
	unless sequences[id][:hom] == adj and sequences[id][:het] == adj
		ratio = sequences[id][:hom].to_f/sequences[id][:het].to_f
	end
	if sequences[id][:hom].to_f > sequences[id][:het].to_f
		sequences[id][:ratio] = ratio
		snpratio[ratio][id] = 1
		ratios << ratio
	end
}

### use maximum ratio value and assingn cutoffs based on it
maximum = ratios.max
minimum = 0
unless filter == 0
	minimum = (maximum * filter) / 100
end

output = File.open("processed_varinfo.txt", "w")
output.puts "fragment\tnumhm\tnumht\tratio"
### select fragments based on cut offs and store ids to an array
selected = []
sequences.keys.each { | id |
	if sequences[id][:ratio] > minimum
		output.puts "#{id}\t#{sequences[id][:hom]}\t#{sequences[id][:het]}\t#{sequences[id][:ratio]}"
		selected << id
	end
}
output.close

rubyR = RinRuby.new(:echo=>false)
rubyR.eval "library(zoo)"
rubyR.eval "source('./peakfind.r')"
rubyR.eval "inputdf = read.delim(file=\"processed_varinfo.txt\", header = TRUE, stringsAsFactors=FALSE)"
selected.permutation.each do | arranged |
	rubyR.arranged = arranged
	rubyR.eval "newdf = inputdf[match(arranged, indf$fragment,]"
	rubyR.eval "row.names(newdf) = NULL"
	rubyR.eval "peaks = argmax(row.names(newdf), newdf$ratio, span=1)"
	numpeak = rubyR.pull "length(peaks$i)"
	if numpeak == 1
		puts "#{auc}\t#{arranged.join("\t")}"
	end
end

rubyR.quit

=begin
new_order = []
n = 1
snpratio.keys.sort.reverse.each { |fraction|
  if fraction >= minimum
  	randomized = snpratio[fraction].keys.shuffle
	randomized.each do |id|
		if n%2 == 0
			new_order.unshift(id)
		else
			new_order.push(id)
		end
		n += 1
	end
  end
}

warn "#{new_order}\n"
=end
