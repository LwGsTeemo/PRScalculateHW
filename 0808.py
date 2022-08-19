import pandas as pd

VCF = open("HWGS10060.vcf", 'r+') # read vcf file
hg38 = pd.read_csv("hg38_pos.csv") # read SMPs excel file
df= pd.DataFrame(hg38)
df = df.dropna(subset=['Chinese OR (95% CI)7'], how='all') # clean the SNP if it doesn't have OR score
chr = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
with open("PRS_score.txt",'w+') as f:
    with open("hg38.txt",'a+') as g:
        PRS=0
        cnt=0
        for line in VCF.readlines():
            if line.startswith('#'): # ignore all header
                continue
            info = line.split("\t")
            if info[0][3:] in chr:
                CHROM = eval(info[0][3:]) # make "chr.1" to "1"
            else:
                CHROM = 0
            POS   = int(info[1]) # get POS
            REF   = info[3] # get REF
            ALT   = info[4] # get ALT
            WGS   = info[-1][0:3] # get WGS eg."0/1" or "1/1"
            for index ,row in df.iterrows():
                if row["Chr."]==CHROM  and row["hg38_pos"]==POS and row["Allele2"].split("/")[0] ==REF and \
                    row["Allele2"].split("/")[1] ==ALT:  # check whether both are the same
                    for i in df.loc[index].tolist(): # write .txt
                        g.write(str(i)+"\t") 
                    g.write("\n")
                    print("find one!")
                    Or = float(df.loc[index]["Chinese OR (95% CI)7"].split(" ")[0])
                    if WGS[0]!='0' and WGS[2]!='0': #calculate PRS score
                        PRS+= 2*Or
                    else:
                        PRS+= Or
                    cnt+=1
                    print("Score:",PRS)
                    break
        print("found ",cnt," .")
        print("PRS score =",PRS)
        f.write("total: "+str(cnt)+"\n")
        f.write("Score: "+str(PRS)+"\n")
    g.close()
f.close()

    


# header style
# CHROM	    POS	     ID	 REF  ALT	QUAL	FILTER	INFO	FORMAT	HWGS10060
# ['chr1', '10904', '.', 'G', 'A', '78.28', '.', 'AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=36.67;QD=25.36;SOR=2.303', 'GT:AD:DP:GQ:PL', '1/1:0,2:2:6:90,6,0\n']