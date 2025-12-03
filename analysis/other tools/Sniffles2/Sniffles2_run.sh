sniffles --input {tumorBam} --snf {tumor.snf} -t {threads} --minsupport 3 --mapq 10 --min-alignment-length 1000 --output-rnames --allow-overwrite --long-ins-length 100000 --reference {reference}  
sniffles --input {controlBam} --snf {control.snf} -t {threads} --minsupport 3 --mapq 10 --min-alignment-length 1000 --output-rnames --allow-overwrite --long-ins-length 100000 --reference {reference} 
sniffles --input {tumor.snf} {control.snf} --vcf {somatic_and_normal_sv.vcf} -t {threads} --minsupport 3 --mapq 10 --min-alignment-length 1000 --output-rnames --allow-overwrite --long-ins-length 100000 --reference {reference}
python extract_somatic_sv.py -i {somatic_and_normal_sv.vcf} -o {somatic_sv.vcf}
