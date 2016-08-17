#!/bin/csh -f

foreach x ( SAM/* )
echo $x
samtools  view -bS $x | samtools sort -o  BAM/`basename -s .sam $x`.bam
end

foreach x ( BAM/*.bam )
echo $x
genomeCoverageBed -strand - -bg  -scale -1 -ibam $x  | awk '{print "chr"$0;  }' | awk '{gsub("chrMT", "chrM", $0); print $0;}' >! WIG/`basename -s .bam $x`-C.bed
genomeCoverageBed -strand + -bg -ibam $x | awk '{print "chr"$0;  }' | awk '{gsub("chrMT", "chrM", $0); print $0;}' >! WIG/`basename -s .bam $x`-W.bed
end


foreach x ( WIG/*.bed )
echo $x
/cs/wetlab/pipeline/tools/bedGraphToBigWig $x sacCer3.sizes /cs/wetlab/Alon/website/asym_nucs_rna/`basename -s .bed $x`.bw
end
