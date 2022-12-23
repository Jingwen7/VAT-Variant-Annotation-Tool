#!~/.conda/envs/lra/bin/python3

"""
Author: Jingwen Ren (jingwenr@usc.edu)
PhD Candidate in Epidemiology
Quantative Biology and Bioinformatics
University of Southern California
"""

import argparse
import os
import sys
import csv
import re
import pysam 
from collections import defaultdict
import vcf
import requests
import json
from requests import HTTPError
import multiprocessing
from multiprocessing.dummy import Pool

# Initialize parser
parser = argparse.ArgumentParser(description="VAT: Variant Annotation Tool")
parser.add_argument("-i", "--Input_file_path", nargs='?', required=True)
parser.add_argument("-o", "--Output_file_path", nargs='?', required=True)
args = parser.parse_args()


"""
Rest client for Ensembl API.
"""
class EnsemblRestClient(object):
    def __init__(self, server="http://grch37.rest.ensembl.org"):
        self.server = server

    def perform_rest_action(self, endpoint, data=None, params=None):
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

        if params:
            endpoint += '?' + urlencode(params)

        try:
            response = requests.post(self.server + endpoint, headers=headers, data=data)
            content = response.text
            if content:
                decode = json.loads(content)

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, data, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return decode
    
    
    def get_vep(self, var_ids):

        data_dict = {"hgvs_notations" : var_ids}
        data = json.dumps(data_dict) 
        response = self.perform_rest_action(endpoint = "/vep/human/hgvs", data = data)
        if response:
            return response
        
        return None


class annotateVariant(object):
    def __init__(self, vcf_path):
        self.vcf_reader = vcf.Reader(open(vcf_path, 'r'))
        self.annotate_res = defaultdict(dict) # cpra_id -> {"TC": , "TR": , }
        self.basicInfo()
    
    """
    get a unique cpra id to the INFO col of each variant 
    """
    def get_cpra_ID(self, record):
        alt_str = ','.join([str(a) for a in record.ALT])
        cpra_id = str(record.CHROM) + ":g." + str(record.POS) + record.REF + ">" + alt_str
        return cpra_id
    
    def basicInfo(self):
        for record in self.vcf_reader:
            if record.FILTER:
                continue 
            
            if record.INFO["TC"] < sum(record.INFO["TR"]):
                continue
                
            cpra_id = self.get_cpra_ID(record)
            self.annotate_res[cpra_id]["CHROM"] = record.CHROM
            self.annotate_res[cpra_id]["POS"] = record.POS
            self.annotate_res[cpra_id]["REF"] = record.REF
            self.annotate_res[cpra_id]["ALT"] = record.ALT
            self.annotate_res[cpra_id]["TC"] = record.INFO["TC"]
            self.annotate_res[cpra_id]["TR"] = sum(record.INFO["TR"])
            self.annotate_res[cpra_id]["GT"] = record.samples[0]['GT']
        return 
    
    """
    decide TYPE of variant
    """
    def svType(self):
        for var in self.annotate_res:

            numAlt = len(self.annotate_res[var]["ALT"])
            types = []
            for i in range(numAlt):
                if len(self.annotate_res[var]["REF"]) == len(self.annotate_res[var]["ALT"][i]) == 1:
                    types.append("snp")
                elif len(self.annotate_res[var]["REF"]) > len(self.annotate_res[var]["ALT"][i]):
                    types.append("del")

                else:
                    types.append("ins")

            self.annotate_res[var]["TYPE"] = types
        
        return 
        
  
    """
    calculate # of reads supporting ALT variants / # of reads supporting REF allele
    assumption:  # of reads supporting REF allele = total reads - # of reads supporting ALT variants
    NOTE: filter out variant where TR of ALT alleles > TC
    """
    def percentReadsSupportALTvsRef(self):
        for var in self.annotate_res:
            try:
                tc = self.annotate_res[var]["TC"]
                tr = self.annotate_res[var]["TR"]
                if tc > tr:
                    self.annotate_res[var]["TF"] = tc - tr
                    self.annotate_res[var]["PRAR"] = tr / (tc - tr) 
                elif tc == tr:
                    self.annotate_res[var]["TF"] = 0
                    self.annotate_res[var]["PRAR"] = 'INF'
            except KeyError:
                print("KeyError")      

    
    def task(self, idx):
        var_list = self.var_chunks[idx]
        client = EnsemblRestClient()
        var_response_list = client.get_vep(var_list)
        return var_response_list
        
            
    """
    Fetch annotation from VEP API
    """
    def fetchVarInFoFromVEP(self, numOfProcessors=min(multiprocessing.cpu_count() - 1, 15)):

        numQueryPerTime = 200
        input_var_list = list(self.annotate_res.keys())
        self.var_chunks = [input_var_list[i:i + numQueryPerTime] for i in range(0, len(input_var_list), numQueryPerTime)]
        chunks = range(len(self.var_chunks))
        numOfProcessors = min(numOfProcessors, len(chunks))

        print("start fetching VEP API with %s processors...."%numOfProcessors)
        pool = Pool(numOfProcessors)
        # try:
        for response in pool.imap_unordered(self.task, chunks):
            
            for i in range(len(response)):
                var_id = response[i]['input']
                self.annotate_res[var_id]["ANNO"] = defaultdict(lambda: 'NA')
                if 'transcript_consequences' not in response[i]:
                    continue 
                numTranscript = len(response[i]['transcript_consequences'])
                av.annotate_res[var_id]['ANNO']['transcript_id'] = ",".join([response[i]['transcript_consequences'][j].get('transcript_id', 'NA') for j in range(numTranscript)])
                av.annotate_res[var_id]['ANNO']['gene_id'] = ",".join([response[i]['transcript_consequences'][j].get('gene_id', 'NA') for j in range(numTranscript)])
                av.annotate_res[var_id]['ANNO']['gene_symbol'] = ",".join([response[i]['transcript_consequences'][j].get('gene_symbol', 'NA') for j in range(numTranscript)])
                av.annotate_res[var_id]['ANNO']["biotype"] = ",".join([response[i]['transcript_consequences'][j].get('biotype', 'NA') for j in range(numTranscript)])
                av.annotate_res[var_id]['ANNO']["impact"] = ",".join([response[i]['transcript_consequences'][j].get('impact', 'NA') for j in range(numTranscript)])
                av.annotate_res[var_id]['ANNO']["gene_symbol"] = ",".join([response[i]['transcript_consequences'][j].get('gene_symbol', 'NA') for j in range(numTranscript)])
                av.annotate_res[var_id]['ANNO']["polyphen_prediction"] = ",".join([response[i]['transcript_consequences'][j].get('polyphen_prediction', 'NA') for j in range(numTranscript)])
                av.annotate_res[var_id]['ANNO']["sift_prediction"] = ",".join([response[i]['transcript_consequences'][j].get('sift_prediction', 'NA') for j in range(numTranscript)])
                av.annotate_res[var_id]['ANNO']["most_severe_consequence"] = response[i].get('most_severe_consequence', 'NA')

                if 'colocated_variants' not in response[i]:
                    continue
                if len(response[i]['colocated_variants']) == 1:
                    av.annotate_res[var_id]['ANNO']["minor_allele"] = response[i]['colocated_variants'][0].get('minor_allele', 'NA')
                    av.annotate_res[var_id]['ANNO']["minor_allele_freq"] = response[i]['colocated_variants'][0].get('minor_allele_freq', 'NA')
                    av.annotate_res[var_id]['ANNO']["rsid"] = response[i]['colocated_variants'][0].get('id', 'NA')  
                    av.annotate_res[var_id]['ANNO']["clin_sig"] = response[i]['colocated_variants'][0].get('clin_sig', 'NA')

                else:
                    av.annotate_res[var_id]['ANNO']["minor_allele"] = response[i]['colocated_variants'][1].get('minor_allele', 'NA')
                    av.annotate_res[var_id]['ANNO']["minor_allele_freq"] = response[i]['colocated_variants'][1].get('minor_allele_freq', 'NA')
                    av.annotate_res[var_id]['ANNO']["rsid"] = response[i]['colocated_variants'][1].get('id', 'NA')
                    av.annotate_res[var_id]['ANNO']["clin_sig"] = response[i]['colocated_variants'][1].get('clin_sig', 'NA')
        # except:
        #     print ("Outer exception caught!")

        pool.close()
        pool.join()
        print("finish fetching VEP API....")
        return 
    
    """
    Write output tsv file
    """
    def out(self, out_file):
        with open(out_file, 'w', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
            writer.writerow(["CHROM", "POS", "REF", "ALT", "GT", "TYPE", "total_read_depth", \
                             "ref_read_depth", "alt_read_depth", "Ratio_supporting_reads_alt_vs_ref", \
                              "rsid", "most_severe_consequence", "minor_allele_freq", \
                             "minor_allele", "clin_sig", "transcript_id", "gene_id", "gene_symbol", "impact", \
                             "biotype", "polyphen_prediction", "sift_prediction"])
            
            for var, info in self.annotate_res.items():
                if 'ANNO' in info and info['ANNO']:
                    writer.writerow([info['CHROM'], info['POS'], info['REF'], ','.join([str(ch) for ch in info["ALT"]]), info['GT'], str(info['TC']), \
                                     str(info['TF']), str(info['TR']), str(info['PRAR']), \
                                     info['ANNO']['rsid'], info['ANNO']['most_severe_consequence'], str(info['ANNO']['minor_allele_freq']), \
                                     info['ANNO']['minor_allele'], str(info['ANNO']['clin_sig']), info['ANNO']['transcript_id'], \
                                     info['ANNO']['gene_id'], info['ANNO']['gene_symbol'], info['ANNO']['impact'], \
                                     info['ANNO']['biotype'], info['ANNO']['polyphen_prediction'], info['ANNO']['sift_prediction'],\
                                    ])
                else:
                    writer.writerow([info['CHROM'], info['POS'], info['REF'], ','.join([str(ch) for ch in info["ALT"]]), info['GT'], str(info['TC']), \
                                     str(info['TF']), str(info['TR']), str(info['PRAR'])] + ['NA'] * 12)
        return 



if  __name__ == "__main__":
    av = annotateVariant(args.Input_file_path)
    av.svType()
    av.percentReadsSupportALTvsRef()
    av.fetchVarInFoFromVEP()
    av.out(args.Output_file_path)
                     
