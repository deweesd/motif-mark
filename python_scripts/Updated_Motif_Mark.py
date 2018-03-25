#!/usr/bin/python

## library imports ## 
import re
import math
import argparse
import random
import cairo


############ argparse commands ############ 
parser = argparse.ArgumentParser(description="Python script displaying binding site motifs. Once motifs are found, using picairo to visualize multiple motifs across multiple sequences (UCSC format).")
parser.add_argument('-f','--fasta_file', help='absolute path to FASTA file of sequences to be searched.', required=True, type=str)
parser.add_argument('-m','--motif_list', help='text file that contains the list of motifs to be visualized, formatted as one motif per line.', required=True, type=str)
args = parser.parse_args()

## Setting the passed variables for argparse ##
fasta_df = parser.fasta_file
Motif_df = parser.motif_list


################################## Cairo_Funtions ########################################


## Iterating through each FASTA seq for "A-Z" letters that are paired ~ reg-ex search ## 
search_exon_seq = re.compile("[A-Z]{1,}")


## Parameters drawn for SVG image in picairo ##
def legend_params_cairo(r, g, b, t, motif_seq):
    context.set_source_rgb(r, g, b)
    context.rectangle(t,7,7,7)
    context.fill()
    context.set_source_rgb(0,0,0)
    context.move_to(t + 10,15)
    context.show_text(motif_seq)
    t = t + 10*(len(motif_seq))
    return(t)


## Given start/end sites in sequence, draws out motifs ##
def motif_line(x_coord,y_coord,start_pos,stop_pos):
    context.move_to(x_coord+start_pos,y_coord)
    context.line_to(x_coord+stop_pos,y_coord)
    context.stroke()


#### Calling all 'DNA IUPAC codes' (referenced from "http://www.bioinformatics.org/sms2/iupac.html") 
## given our motif_list ~ will be used in regex code below: ####

def adjusting_motifs(motif_seq):
    motif_seq = motif_seq.upper() #converting to all upper_case letters
    motif_seq = motif_seq.replace('N','.') # any base
    
    #3_base_total_seq ~ Replacing all single_IUPAC code(i.e., B, D, etc.) with [Nucl_base_seqs]
    motif_seq = motif_seq.replace('B','[CGT]') 
    motif_seq = motif_seq.replace('D','[AGT]')
    motif_seq = motif_seq.replace('H','[ACT]')
    motif_seq = motif_seq.replace('V','[ACG]')
    
    #2_base_total_seq ~ Replacing all single_IUPAC code(i.e., Y, R, S, etc.) with [Nuc_base_seqs]
    motif_seq = motif_seq.replace('Y','[CT]')
    motif_seq = motif_seq.replace('R','[AU]')
    motif_seq = motif_seq.replace('S','[GC]')
    motif_seq = motif_seq.replace('W','[AT]')
    motif_seq = motif_seq.replace('K','[GT]')
    motif_seq = motif_seq.replace('M','[AC]')
    return motif_seq



#################### Begin iterating through both FASTA file and motif list ############## 
#################### to establish motifs found in said FASTA file ########################

## Tally total number of lines in the file ##
fasta_df = open('fasta_sample.txt')
for seq_line, value in enumerate(fasta_df, 1):
    print(seq_line, value)
fasta_df.close()

## Total number of lines (including header) ~ 39
total_line_cnt = 0
fasta_df = open('sequences.txt')
for sequence_lines in fasta_df:
    total_line_cnt+=1 # increment by 1
    #print(sequence_lines)
fasta_df.close()

## Make an iterable list from our fasta_df ##
fasta_df = open('sequences.txt')
seq = ""
fasta_list = []
nmbr_cnt=0
for sequence_lines in fasta_df:
    sequence_lines = sequence_lines.rstrip() # chars in each line have been stripped from the end of the string 
    nmbr_cnt+=1
    
    # Take only first header from fasta file and add to fasta_list
    if nmbr_cnt == 1:
        fasta_list.append(sequence_lines)
        seq = ""
        #print(fasta_list)
        
    
    # Look at all 39 lines and add last seq (line 39) and append seq to fasta list
    elif nmbr_cnt == total_line_cnt:
        seq = seq + sequence_lines
        fasta_list.append(seq)
        #print(fasta_list)
    
    # Looking at only headers in fasta file and adding to list
    elif sequence_lines.startswith(">") : # calling header line 
        #append seq and each line to fasta list given header_line == TRUE
        fasta_list.append(seq)
        fasta_list.append(sequence_lines)
        seq = ""
        #print(fasta_list)
    
    # Looking at only sequnces in fasta, ignoring headers 
    elif not sequence_lines.startswith(">"): # ONLY seq lines 
        seq = seq + sequence_lines
        #print(fasta_list)

fasta_df.close()


## Open our fasta file to count the number of genes to find the longest present and account for plot Y axis
entry_len = 0
nmbr_lines = 0
for sequence_lines in fasta_list:
    
    if not sequence_lines.startswith(">"):
        
        nmbr_lines+=1 # Increment by 1
        
        # If len of lines is greater than entry_lent, set entry_len equal to number of lines
        if len(sequence_lines) > entry_len:
            entry_len = len(sequence_lines)
            #print(nmbr_lines)
            


## Take motif_seqs and apply to function while appending to open list ## 
Motif_df = open('Motif_list.txt')
motif_list = []
for motif_seq in Motif_df:
    
    motif_seq = motif_seq.rstrip() 
    #print(motif_seq)
    
    # append converted motifs from 'adjusting_motifs' function above given motif_seq in Motif_df 
    motif_list.append(adjusting_motifs(motif_seq))
    #print(motif_list)

Motif_df.close()


## Generate lists given color selection that corresponds to each motif found in FASTA file ##
green = []
blue = []
red = []

## for each motif seq in our motif_list, randomly assign int_values to said lists
for motif_seq in motif_list:
    green.append(random.randint(0,100)/100)
    #print(green)
    
    red.append(random.randint(0,100)/100)
    
    blue.append(random.randint(0,100)/100)



## Surface coordinates ##
x_coord = 30 # Set x_coord_start value 
y_coord = 0 # Set y_coord_start value

## set surface to SVG_output file while also setting surface given calculated values from x_cord/entery_len/len of seq
surface = cairo.SVGSurface("test.svg", entry_len+(x_coord*2), (nmbr_lines*50)+50)
context = cairo.Context(surface)
context.set_font_size(5) # Font_size can change given size of FASTA file/output 


## setting start position within pycairo setup ##
surface = cairo.SVGSurface("Leslie_plot_2.svg", 1350, 1000)
context = cairo.Context(surface)
context.set_line_width(0.5)
start = [300,100]

## Setting motif_color ## 
t = 30

# call rgb given the motif_list
for motif_seq, r, g, b in zip(motif_list, red, green, blue):
    
    if t == 30:
        t = legend_params_cairo(r,g,b,t,motif_seq)
        #print(n)
    else:
        t = legend_params_cairo(r,g,b,t,motif_seq)
        #print(n)


################### Begin to set the intron/exon positions given sequence and ############ 
################### plot each motif on said sequences using Picairo ######################


for sequence_lines in fasta_list:
    ## look at headers in sequence per chrom_pos ##
    if sequence_lines.startswith(">"): 
        
        # set coord to y_value + 60 for range set
        y_coord = y_coord + 50 
        
        # Set coord of motif_list to 0
        context.set_source_rgb(0,0,0)
        
        context.move_to(x_coord, y_coord-15)
        
        # plot text of seq_line headers
        context.show_text(sequence_lines)
    
    ## now look at sequences and plot intron/exons ##
    if not sequence_lines.startswith(">"): 
        # width is constant (0.5) ~ can change if desirable
        context.set_line_width(0.5)
        context.move_to(x_coord, y_coord)
        context.line_to(x_coord+len(sequence_lines),y_coord)
        context.stroke()
        
        # line width = 8
        context.set_line_width(8)
        
        # Set and reset our exon start and end variables for each sequence
        start_pos_exon = []
        end_pos_exon = []
        
        # Save the location of every uppercase character in the sequence
        for nuc_bases in search_exon_seq.finditer(sequence_lines): 
            
            # 0_based position factors in start_pos plus 1
            start_pos_exon.append(nuc_bases.start()+1) 
            # end_based position given nuc_bases @ end of seq
            end_pos_exon.append(nuc_bases.end())
        
        # Iterate through our start and stop sites to draw our exons
        for start_pos, stop_pos in zip(start_pos_exon, end_pos_exon): 
            motif_line(x_coord,y_coord,start_pos,stop_pos)
            #print(start_pos)
        
        # Setting seq to "UPPER" ~ looking for uppercass_bases
        sequence_lines = sequence_lines.upper()
        
        # Go through motif_list while factoring the color_scheme set for corresponding motif
        for motif_seq, r, g, b in zip(motif_list, red, green, blue):
            
            
            # Empty list for motifs start/end locations
            start_pos_motif = []
            end_pos_motif = []
            
            motifs = re.compile(motif_seq,) 
        
            # Return an iterator over all non-overlapping matches for the RE pattern in our sequences
            for nuc_bases in motifs.finditer(sequence_lines):
                # 0_based position factors in start_pos plus 1
                start_pos_motif.append(nuc_bases.start()+1) 
                end_pos_motif.append(nuc_bases.end())
        
            # Draw motifs given positions found
            for start_pos, stop_pos in zip(start_pos_motif, end_pos_motif):         
                context.set_source_rgb(r, g, b)
                motif_line(x_cord,y_cord,start_pos,stop_pos)

#write out drawing
surface.finish()
        






