#/bin/bash
main=""
scripts=$(dirname "`realpath $0`")

#####
mkdir ILS_introgression

#QuIBL
mkdir -p QuIBL/inout
lengthtree="$main/phylogeny/CDS_align_merge.fa.raxml.bestTree"
ln -s $Iqtree QuIBL/
ln -s $lengthtree QuIBL/snake.branchlength.nwk
echo -e "[Input]\n../iqtree.all.tree\nnumdistributions: 2\nlikelihoodthresh: 0.01\nnumsteps: 50\ngradascentscalar: 0.5\ntotaloutgroup: Adia\nmultiproc: True\nmaxcores: 10\n\n[Output]\nOutputPath: ./Out.csv" >QuIBL/inout/sampleInputFile.txt
cd QuIBL/inout
python QuIBL.py sampleInputFile.txt 1>1.log 2>2.log #https://github.com/miriammiyagi/QuIBL
cd .. && mkdir result && cd result
python $scripts/Format_QuIBL_final.py ../inout ../snake.branchlength.nwk Lnig
cd ../../

##ABBA-BABA test 
mkdir dfoil 
dfoilpath="" #https://github.com/jbpease/dfoil
MergCDS_alignemt="$main/phylogeny/CDS_align_merge.fa"
echo "node\truns\tp1\tp2\tp3\tpo\tD-stats\tP-value\tmodel" >$main/ILS_introgression/dfoil/stats.txt 
cat $main/introgression_run.list |while read nodeid runs p1 p2 p3 po 
do
        runout="$main/ILS_introgression/dfoil/$nodeid/$runs"
        mkdir -p $runout
        cd  $runout
        for i in p1 p2 p3 po ; do samtools faidx $MergCDS_alignemt $i ;done >$runout/cds.gap.fa #https://github.com/samtools/samtools
        trimal -in $runout/cds.gap.fa -noallgaps -out $runout/cds.ungap.fa -fasta #https://github.com/inab/trimal
        python $dfoilpath/fasta2dfoil.py $runout/cds.ungap.fa -o $runout/cds.gap.fa.count --names $p1,$p2,$p3,$po 1>>1.log 2>>2.log
        python $dfoilpath/dfoil.py --infile $runout/cds.gap.fa.count --out $runout/cds.gap.fa.count.dfoil --mode dstat --plot $runout/plot.pdf 1>>1.log 2>>2.log
        reval=`grep -v "#" cds.gap.fa.count.dfoil|awk '{print $10"\t"$11"\t"$12}'` >>$main/ILS_introgression/dfoil/stats.txt
done
