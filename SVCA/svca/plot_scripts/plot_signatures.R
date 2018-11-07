#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 2){
  print("USAGE: please give the SVCA output directory as a first argument, output plot directory as a second argument")
}
source('visu_signatures_functions.R')

working_dir = args[1]
plot_dir = args[2]
plot_sig_all(working_dir, plot_dir)
