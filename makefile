PredictionUsingChIPSeq.txt:predictionBasedChIPseq.R
	Rscript predictionBasedChIPseq.R

%.fa.smat_8mer.gz: %.fa fasta2matrix_sparse.py
	python fasta2matrix_sparse.py -mismatch 3 -upto -revcomp -normalize frequency 8  $<

%.fa.5merPair.smat.gz: %.fa fasta2sigKmerPair_batch.py
	python fasta2sigKmerPair_batch.py -mismatch 0 -upto -revcomp -normalize frequency 5 $<

%.fa.5merPair.vistasel.smat.gz:%.fa data/vista.fa.5merPair.smat.gz.name fasta2sigKmerPair_batch_sel.py
	python fasta2sigKmerPair_batch_sel.py -upto -revcomp -normalize frequency  -kmerfile data/vista.fa.5merPair.smat.gz.name 5 $<

%.fa.name:%.fa
	grep ">" $< |sed 's/>//' > $@
#PredictionUsingKmer.txt:predictionBasedKmer.R data/LBNLtest.mm10.fa.smat_8mer.gz  data/LBNLtest.mm10.fa.5merPair.vistasel.smat.gz data/LBNLtest.mm10.fa.name
#	Rscript $<

GenomeWide5kPredictionUsingKmer.txt:GWpredictionBasedKmer.R data/GW5k.mm10.fa.smat_8mer.gz  data/GW5k.mm10.fa.5merPair.vistasel.smat.gz data/GW5k.mm10.fa.name
	Rscript GWpredictionBasedKmer.R

README.pdf:README.md
	pandoc -V geometry:"width=25cm , margin=0.1in"  -s  README.md -o README.pdf

data/GW5k.mm10.fa:GenomeWide5kPredictionUsingChIPSeq.txt
	tail -n+2 GenomeWide5kPredictionUsingChIPSeq.txt|tr '|' '\t'|awk '{m=int(($$2+$$3)/2);OFS="\t";print $$1,m-500,m+500,$$1"|"$$2"|"$$3}'|fastaFromBed -fi ~/archive/genome/mm10.fa -bed stdin -name -fo $@ 




submission/set1_file1.bed5:PredictionUsingChIPSeq.txt data/heart.bed 
	tail -n+2 PredictionUsingChIPSeq.txt|cut -f 1,2|tr '|' '\t' |windowBed -u -w 1 -a stdin -b $(word 2, $^)> $@

submission/set1_file2.bed5:PredictionUsingChIPSeq.txt  data/fb.bed
	tail -n+2 PredictionUsingChIPSeq.txt|cut -f 1,3|tr '|' '\t' |windowBed -u -w 1 -a stdin -b $(word 2, $^)> $@	

submission/set1_file3.bed4:GenomeWide5kPredictionUsingChIPSeq.txt
	tail -n+2 $< |awk '$$4>0.8{OFS="\t";print $$1,$$4}' |tr '|' '\t' > $@

submission/set1_file4.bed4:GenomeWide5kPredictionUsingChIPSeq.txt
	tail -n+2 $< |awk '$$3>0.8{OFS="\t";print $$1,$$3}' |tr '|' '\t' > $@



submission/set2_file1.bed5:PredictionUsingKmer.txt data/heart.bed
	tail -n+2 PredictionUsingChIPSeq.txt|cut -f 1,2|tr '|' '\t' |windowBed -u -w 1 -a stdin -b $(word 2, $^)> $@

submission/set2_file2.bed5:PredictionUsingKmer.txt  data/fb.bed
	tail -n+2 PredictionUsingChIPSeq.txt|cut -f 1,3|tr '|' '\t' |windowBed -u -w 1 -a stdin -b $(word 2, $^)> $@    

submission/set2_file3.bed4:GenomeWide5kPredictionUsingKmer.txt
	tail -n+2 $< |awk '$$4>0.5{OFS="\t";print $$1,$$4}' |tr '|' '\t' > $@

submission/set2_file4.bed4:GenomeWide5kPredictionUsingKmer.txt
	tail -n+2 $< |awk '$$3>0.5{OFS="\t";print $$1,$$3}' |tr '|' '\t' > $@


submission.%:submission/set%_file1.bed5 submission/set%_file2.bed5 submission/set%_file3.bed4 submission/set%_file4.bed4

submission:submission.1 submission.2
