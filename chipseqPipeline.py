__author__ = 'lucassilva'
# -*- coding: utf-8 -*-
import sys
import os
from glob import glob
import re


def create_config_file():
    config='[configurações_gerais]\nprocessadores=\nindex_bowtie=\nchr_size_file=\n\n<marca1>\nnome_marca=\n' \
           'sra_file=\n\n<marca1>\nnome_marca=\nsra_file=\n\n\n\n#adicione novas entradas seguindo o mesmo esquema\n' \
           '#   <marcaX>\n#   nome_marca=\n# sra_file=\n#\n'

    salvar_config = open('chip_seq_pipeline.cfg', 'w')
    salvar_config.write(config)
    salvar_config.close()

    print config

def sraToFastq(sra, diretorio):

    #os.system('pwd')
    #exit()
    if not diretorio in glob("*"):

        os.system('mkdir {}'.format(diretorio))
        os.chdir('{}'.format(diretorio))

        if re.search('fa$|fastq$', sra.split('.')[-1]):
            sys.stderr.write("Criando link simbolico do arquivo fa/fastq..." + "\n")
            os.system('ln -T {sra} {sra_link}'.format(sra=sra, sra_link=sra.split('/')[-1]))
        elif re.search('sra$', sra.split('.')[-1]):
            sys.stderr.write('Descompactando arquivo SRA...\n')
            os.system('fastq-dump {}'.format(sra))
        else:
            sys.stderr.write("extension error..." + "\n")
    else:
        sys.stderr.write('Diretorio já existente....\nEntrando no Diretório\n')
        print diretorio
        #os.system('pwd')
        os.chdir('{}'.format(diretorio))
        #os.system('pwd')
        if re.search('fa$|fastq$', sra.split('.')[-1]):
            sys.stderr.write("Criando link simbolico do arquivo fa/fastq..." + "\n")
            os.system('ln -T {sra} {sra_link}'.format(sra=sra, sra_link=sra.split('/')[-1]))
        elif re.search('sra$', sra.split('.')[-1]):
            sys.stderr.write('Descompactando arquivo SRA...\n')
            os.system('fastq-dump {}'.format(sra))
        else:
            sys.stderr.write("extension error..." + "\n")

def controle_qualidade(sra):
    sys.stderr.write('Tratando com filtros de qualidade o arquivo fastq...\n')
    print sra.split('.')[-1]
    if sra.split('.')[-1].lower() == 'sra':
        sys.stderr.write("Preparando FASTQ para filtragem..." + "\n")
        nome_arquivo_sra = sra.split('/')[-1].replace(".sra", ".fastq")

    if sra.split('.')[-1].lower() == 'fa':
        nome_arquivo_sra = sra.split('/')[-1]
        #nome_fasta_filtrado= 'filtro_tamanho_trimado_{}'.format(nome_arquivo_sra)
        nome_fasta_filtrado= nome_arquivo_sra
        sys.stderr.write("Arquivo .fa encontrado. Pulando etapa de filtragem..." + "\n")
        return nome_fasta_filtrado
    if sra.split('.')[-1].lower() == 'fastq':
        nome_arquivo_sra = sra.split('/')[-1]
        #nome_fasta_filtrado= 'filtro_tamanho_trimado_{}'.format(nome_arquivo_sra)
        nome_fasta_filtrado= nome_arquivo_sra
        sys.stderr.write("Preparando FASTQ para filtragem..." + "\n")


    #filtrar bases com qualidade baixa
    os.system('fastqutils filter -qual 30 3 {} > trimado_{}'.format(nome_arquivo_sra, nome_arquivo_sra))
    # filtrar bases que ficaram menor que 30
    os.system('fastqutils filter -size 30 trimado_{} > filtro_tamanho_trimado_{}'.format(nome_arquivo_sra, nome_arquivo_sra))
    nome_fasta_filtrado= 'filtro_tamanho_trimado_{}'.format(nome_arquivo_sra)
    return nome_fasta_filtrado

def mapeamento(index, arquivo_fasta_filtrado, processadores):
    sys.stderr.write('Mapeando arquivos contra o genoma de referencia...\n')
    if arquivo_fasta_filtrado.split('.')[-1].lower() == 'fastq':
        out_sam= arquivo_fasta_filtrado.replace('.fastq', '.sam')
        comando='bowtie2 -x {index_file} -p {processadores} {arquivo_fasta_filtrado} > {saida_sam}'.format(index_file = index,
                                                                                              processadores = processadores,
                                                                                      arquivo_fasta_filtrado =arquivo_fasta_filtrado,
                                                                                              saida_sam=out_sam)

    if arquivo_fasta_filtrado.split('.')[-1].lower() == 'fa':
        out_sam= arquivo_fasta_filtrado.replace('.fa', '.sam')
        comando='ls ;bowtie2 -x {index_file} -p {processadores} -f {arquivo_fasta_filtrado} > {saida_sam}'.format(index_file = index,
                                                                                              processadores = processadores,
                                                                                      arquivo_fasta_filtrado =arquivo_fasta_filtrado,
                                                                                              saida_sam=out_sam)


    os.system(comando)
    return out_sam

def criando_diretorio_tag(sam_file):
    sys.stderr.write('Criando diretorio tag para utilizar o homer...\n')
    #criando diretorio com tags
    os.system('makeTagDirectory tag_dir {}'.format(sam_file))

def find_peaks():
    sys.stderr.write('Encontrando picos nos dados de chip-seq...\n')
    '''UCSC
    Find Peaks - Filtrando baseado no local signal onde o -L indica ter densidade de pico maior que 4x do que a região de 10kb
    o -center irá centralizar os reads com overlap, aplicando-se isso para marcas shark'''

    os.system('findPeaks tag_dir -center -L -o peakfile > peakfile')

def convertendo_para_outros_formatos(saida_name, chr_size):
    sys.stderr.write('Convertendo picos para bed, bedGraph e bigWig...\n')
    sys.stderr.write('....BED....\n')
    os.system('pos2bed.pl peakfile > {}.bed'.format(saida_name))
    sys.stderr.write('....BEDGRAPH...\n')
    #Gerar o bedGraph
    os.system('makeUCSCfile  tag_dir -o {}'.format(saida_name))
    sys.stderr.write('....BigWig...\n')
    #Gerar o bigWig
    print ('makeUCSCfile tag_dir -o {}.bigWig -bigWig {} -fsize 1e20 > trackinfo.txt'.format(saida_name, chr_size))
    os.system('makeUCSCfile tag_dir -o {}.bigWig -bigWig {} -fsize 1e20 > trackinfo.txt'.format(saida_name, chr_size))

def removendo_arquivos_temp_criados():
    '''remove arquivos criados durante o processo'''
    sys.stderr.write('Removendo todos os arquivos *fastq/*fa | *sam e o diretorio tag_dir do diretorio\n')
    os.system('rm *fastq ; rm *sam ; rm *fastq ; rm *fa; rm -r tag_dir')

def main():

    run_dir = os.getcwd()
    if '-create_cfg'in sys.argv:
        create_config_file()
        sys.stderr.write('arquivo <chip_seq_pipeline.cfg> criado no config criado no diretorio atual\n')
        exit()
    if '-config' in sys.argv:
        try:
            file_index = sys.argv[sys.argv.index('-config')+1]
            arquivo = open(file_index, 'r').read()
            #parse opções globais
            #for x in re.split('\[\S*\]\n',arquivo):
            filtrando=re.sub('^\s+| +', '', arquivo)
            arquivo = re.sub('#\S+\n','', filtrando)
            #print arquivo
            config_geral = re.split('\[configurações_gerais\]\n', arquivo)[1].split('<')[0]
            parse_geral = [valor for valor in config_geral.split('\n') if valor and '#' not in valor]
            #print parse_geral
            pc, ib, cs = parse_geral
            if "" in [elemento.split('=')[1] for elemento in [pc, ib, cs]]:
                sys.stderr.write("campos vazios ou corrompidos encotrados no campo de CONFIGURAÇÕES GERAIS  arquivo de configuração.\n")
                exit()

            amostras=''

            if len(re.split('\<\S+>\n', arquivo)[1:]) == 1:

                amostras=re.split('\<\S+>\n', arquivo)[0]
            else:
                amostras=re.split('\<\S+>\n', arquivo)[1:]

            print amostras
            for amostra_x in amostras:

                os.chdir(run_dir)

                nome_marca, sra_nome = [r_amostra for r_amostra in  amostra_x.split('\n') if r_amostra and '#' not in r_amostra]

                print nome_marca, sra_nome
                temp_fastq_filtrado=''
                diretorio = nome_marca.split('=')[1]
                nome_diretorio = nome_marca.split('=')[1].split('/')[-1]
                sra_file = sra_nome.split('=')[1]
                bowtie_index = ib.split('=')[1]
                processadores = pc.split('=')[1]
                chr_size_file = cs.split('=')[1]

                if '' in [sra_file, diretorio] :
                    sys.stderr.write("campos vazios ou corrompidos encotrados no campo das AMOSTRAS(sra_name ou nome_marca)  arquivo de configuração.\n")
                sys.stderr.write('Rodando amostras para arquivo {}\n'.format(diretorio))

                #print diretorio, sra_file, bowtie_index
                ##############pipeline#####################
                sraToFastq(sra_file, nome_diretorio)
                fasta_gerado = controle_qualidade(sra_file)
                arquivo_sam = mapeamento(bowtie_index, fasta_gerado, processadores)
                criando_diretorio_tag(arquivo_sam)
                find_peaks()
                convertendo_para_outros_formatos(nome_diretorio, chr_size_file)

                if '-save_temp'in sys.argv:
                    pass
                else:
                    try:
                        removendo_arquivos_temp_criados()
                    except:
                        pass

        except:
            sys.stderr.write('Erro:\n________\nA sintaxe correta é -config <caminho do arquivo de configuração>\nCaso não exista'
                             ' um arquivo de configuração (ou esse arquivo seja danificado) utilize \"python chipseqPipeline.py -create_cfg  \" para criar um novo. \n')
            exit()
    #caso não exista uma arquivo de configuração
    else:

        if  len(sys.argv) != 11:
            sys.stderr.write('----ERRO---!\n Entre com os seguintes parametros\n-d  nome do diretorio a ser criado \n-sra nome do arquivo sra utilizado na analise (ou arquivo fastq/fa)'
                             '\n-b index do bowtie\n-p numero de processadores utilizados no alinhamento \n-cs arquivo contendo o tamanho dos cromossomos \n'
                             'Use -config <arquivo de configuração> para rodar em lote (utilize -create_cfg para criar um arquivo novo de configuração)\n')
            exit()


        print "começando programa"
        temp_fastq_filtrado=''
        diretorio = sys.argv[sys.argv.index('-d')+1]
        nome_diretorio = diretorio.split('/')[-1]
        sra_file = sys.argv[sys.argv.index('-sra')+1]
        bowtie_index = sys.argv[sys.argv.index('-b')+1]
        processadores = sys.argv[sys.argv.index('-p')+1]
        chr_size_file = sys.argv[sys.argv.index('-cs')+1]
        print diretorio, sra_file, bowtie_index
        ##############pipeline#####################
        sraToFastq(sra_file, nome_diretorio)
        fasta_gerado = controle_qualidade(sra_file)
        arquivo_sam = mapeamento(bowtie_index, fasta_gerado, processadores)
        criando_diretorio_tag(arquivo_sam)
        find_peaks()
        convertendo_para_outros_formatos(nome_diretorio, chr_size_file)

        if '-save_temp'in sys.argv:
            pass
        else:
            try:
                removendo_arquivos_temp_criados()
            except:
                pass


if __name__ == '__main__':
    sys.exit(main())

