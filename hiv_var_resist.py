#!/usr/bin/env python
"""Variant call HIV deep sequenced populations, assigning drug resistance mutations.

- Start with from BAM inputs files from bcbio following disambiguation
  to remove human from HIV reads.
- Runs hivmmer (https://github.com/kantorlab/hivmmer) to call variants
  on POL genes, producing AAVF output.
- Annotates with drug resistance mutations from Stanford HIVdb. Uses
  PyAAVF parser (https://github.com/winhiv/PyAAVF) for AAVF inputs,
  ASIpython (https://github.com/phac-nml/asipython) for parsing and evaluating
  XML from Stanford HIVdb. aavf-tools (https://github.com/winhiv/aavf-tools)
  is a useful code base for learning how to use these.

Usage:
    hiv_var_resist.py
      <input CSV with flowcell and sample names>
      <flowcell directory containing bcbio outputs>
      <HMM reference for gene of interest, from hivmmer GitHub repo>
      <HIVdb XML file, from aavf-tools GitHub repo>
"""
import collections
import csv
import functools
import operator
import os
import subprocess
import sys

from PyAAVF import parser, model
from Asi.XML.XmlAsiTransformer import XmlAsiTransformer
from Asi.Grammar.StringMutationComparator import StringMutationComparator

def main(sample_csv, flowcell_dir, reference_hmm, hivdb_xml):
    summarize_hivdb(hivdb_xml)
    calls = []
    with open(sample_csv) as in_handle:
        for flowcell, sample in (l.strip().split(",") for l in in_handle if l.strip()):
            print(flowcell, sample)
            bam_file = os.path.join(flowcell_dir, flowcell, flowcell, "final", sample,
                                    "%s-ready.bam" % sample)
            work_dir = os.path.join("work", sample)
            if not os.path.exists(work_dir):
                os.makedirs(work_dir)
            calls.append((sample, run_sample(sample, bam_file, reference_hmm, hivdb_xml, work_dir)))
    summarize_sample_calls(sample_csv, calls)

def run_sample(sample, bam_file, reference_hmm, hivdb_xml, work_dir):
    r1 = os.path.join(work_dir, "%s_R1.fq" % (sample))
    r2 = os.path.join(work_dir, "%s_R2.fq" % (sample))
    cmd = ["bamtofastq", "filename=%s" % bam_file, "F=%s" % r1, "F2=%s" % r2,
           "S=/dev/null", "O=/dev/null", "O2=/dev/null"]
    if not os.path.exists(r1):
        subprocess.check_call(cmd)
    aavf_file = os.path.join(work_dir, "%s-hivmmer.hmmsearch2.aavf" % sample)
    if not os.path.exists(aavf_file):
        cmd = ["hivmmer", "--id", os.path.join(work_dir, "%s-hivmmer" % sample),
              "--fq1", r1, "--fq2", r2, "--ref", reference_hmm]
        subprocess.check_call(cmd)
    aavf_res_file = aavf_file.replace(".aavf", ".annotated.aavf")
    if not os.path.exists(aavf_res_file):
        annotate_aavf(aavf_file, hivdb_xml, aavf_res_file)
        # cmd = ["aavfresistance", "-a", aavf_file, "-x", hivdb_xml, "-o", aavf_res_file]
        # subprocess.check_call(cmd)
    print(aavf_res_file)
    return aavf_res_file

def evaluate_resistance(recs, genes):
    """Evaluate drug resistance for the given set of variant records.
    """
    mutations = collections.defaultdict(list)
    for rec in recs:
        mutations[rec.GENE].append("%s%s%s" % (rec.REF, rec.POS, rec.ALT[0]))

    rmuts = collections.defaultdict(lambda: collections.defaultdict(list))
    for gene in mutations:
        eval_gene = genes[gene].evaluate(mutations[gene], StringMutationComparator(True))
        for drug_class in eval_gene.get_evaluated_drug_classes():
            for drug in drug_class.get_evaluated_drugs():
                is_resistant = False
                for c in drug.get_evaluated_conditions():
                    for d in c.get_definitions():
                        if d.sir == "R":
                            is_resistant = True
                            break
                if is_resistant:
                    for mut in drug.get_scored_mutations():
                        key = (gene, mut)
                        rmuts[key][drug_class.get_drug_class().name].append(drug.get_drug().name)
    return rmuts

def annotate_aavf(in_file, hivdb_file, out_file):
    """Annotate an AAVF input file with drug resistance changes.
    """
    res_genes = XmlAsiTransformer(True).transform(open(hivdb_file, "r"))

    reader = parser.Reader(in_file)
    aavf_obj = reader.read_records()
    aavf_obj.infos["CAT"] = model.Info("CAT", ".", "String", "Drug resistance category", None, None)
    aavf_obj.infos["DRUG"] = model.Info("DRUG", ".", "String", "Drug reistances", None, None)
    with open(out_file, "w") as out_handle:
        writer = parser.Writer(out_handle, aavf_obj)
        for rec in aavf_obj:
            if rec.POS not in rec.ALT:
                rmuts = evaluate_resistance([rec], res_genes)
                k = (rec.GENE, "%s%s%s" % (rec.REF, rec.POS, rec.ALT[0]))
                if rmuts.get(k):
                    cats = []
                    drugs = []
                    for cat in rmuts.get(k).keys():
                        cats.append(cat)
                        drugs.extend(rmuts[k][cat])
                    rec.INFO["CAT"] = cats
                    rec.INFO["DRUG"] = drugs
                    writer.write_record(rec)

def summarize_sample_calls(in_file, aavf_files):
    """Provide a single tab delimited output of samples from AAVF inputs.
    """
    out_file = os.path.join(os.getcwd(), "%s-calls.csv" % os.path.splitext(os.path.basename(in_file))[0])
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["sample", "gene", "pos", "ref", "alt", "altfreq", "altcount", "coverage",
                         "filter", "drug", "drugcat"])
        for sample, aavf_file in aavf_files:
            reader = parser.Reader(aavf_file)
            for rec in reader.read_records():
                writer.writerow([sample, rec.GENE, rec.POS, rec.REF, ";".join(rec.ALT),
                                 rec.ALT_FREQ, rec.COVERAGE,
                                 int(functools.reduce(operator.add, rec.INFO.get("ACC", []))),
                                 ";".join(rec.FILTER) if rec.FILTER else "pass",
                                 ";".join(rec.INFO.get("DRUG", [])), ";".join(rec.INFO.get("CAT", []))])
    return out_file

def summarize_hivdb(hivdb_file):
    """Output genes and conditions covered in HIVdb XML output.
    """
    res_genes = XmlAsiTransformer(True).transform(open(hivdb_file, "r"))
    with open("hivdb_summary.txt", "w") as out_handle:
        for gene in res_genes.values():
            for drug_class in gene.get_drug_classes():
                for drug in drug_class.get_drugs():
                    for rule in drug.get_drug_rules():
                        out_handle.write("%s %s %s %s\n" % (gene.name, drug_class.name,
                                                            drug.name, rule.condition.statement))

if __name__ == "__main__":
    main(*sys.argv[1:])
