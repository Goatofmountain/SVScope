#!/usr/bin/env python3
import argparse


def sniffles2_somatic_sv_searching(multisample_vcf, somatic_vcf):
    with open(multisample_vcf, 'r') as mvf, open(somatic_vcf, 'w') as svf:
        for line in mvf:
            if line.startswith('#'):  
                svf.write(line)
            else:
                content = line.strip().split('\t')  
                info = content[7].split(';')
                SUPP_VEC = info[-1].split("=")[1]
                if SUPP_VEC == '10':
                    svf.write(line)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('-i', '--input', required=True, help='multi-sample VCF file')
    parser.add_argument('-o', '--output', required=True, help='output somatic VCF file')
    args = parser.parse_args()
    sniffles2_somatic_sv_searching(args.input, args.output)
