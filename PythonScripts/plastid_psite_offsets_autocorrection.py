#! /usr/bin/python
#-*- coding: utf-8 -*-

import argparse
import sys
import os
import math

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", required = True, metavar = "original p site offset file", help = "original p site offset file to correct")
	parser.add_argument("-m", "--metageneprofiles", required = True, metavar = "metagene profiles", help = "the metagene profiles")
	parser.add_argument("-o", "--output", required = True, metavar = "corrected p site offset file")
	args = parser.parse_args()
	
	output_poffset_corrected = args.output
	
	if os.path.exists(args.input):
		input_poffset = args.input
	else:
		print "No such file : " + args.input
		sys.exit(1)
	
	if os.path.exists(args.metageneprofiles):
		input_metagene_profiles = args.metageneprofiles
	else:
		print "No such file : " + args.metageneprofiles
		sys.exit(1)
	
	i = open(input_poffset, "r")
	m = open(input_metagene_profiles, "r")
	o = open(output_poffset_corrected, "w")
	
	dict_original_poffset = {}
	dict_final_poffset = {}
	dict_metagene_profiles = {}
	
	for line in m:
		if not line.startswith("#"):
			if line.startswith("x"):
				l = line.strip().split()
				lgt = l[1:]
				lgt = [c.rstrip("-mers") for c in lgt]
				
			else:
				l = line.strip().split()
				offset = int(l[0])
				l_figures = l[1:]
			
				for f in range(len(l_figures)):
					if offset <= 0:
						if lgt[f] not in dict_metagene_profiles:
							dict_metagene_profiles[lgt[f]] = []
						
						dict_metagene_profiles[lgt[f]].append((offset, float(l_figures[f])))
	
	#print dict_metagene_profiles
	
	for line in i:
		if not line.startswith("#") and not line.startswith("length") and not line.startswith("default"):
			l = line.strip().split()
			dict_original_poffset[l[0]] = int(l[1])
	
	dict_final_poffset = dict_original_poffset
	
	#print list(reversed(sorted(dict_original_poffset.items())))
	
	l_poffset = list(reversed(sorted(dict_original_poffset.items())))
	
	for x in range(len(l_poffset)):
		lgt1 = l_poffset[x]
		
		if x < len(l_poffset) - 1:
			lgt2 = l_poffset[x + 1]
		
		#if lgt1 == lgt2:
		#	lgt = lgt1
		#	lgt2 = None
		#	lgt1 = None
			
		if x > 0:
			lgt_previous = l_poffset[x - 1]
			print lgt_previous
			diff_poffset_current_previous = math.fabs(int(lgt1[1]) - int(lgt_previous[1]))
		else:
			diff_poffset_current_previous = None
			
		if x == len(l_poffset) - 1:
			lgt2 = lgt_previous
		
		#if lgt1 != None and lgt2 != None:
		#print x
		print lgt1
		print lgt2
		
		if math.fabs(int(lgt1[1]) - int(lgt2[1])) > 1:
			print "WEIRD"
			lowest_poffset = min(int(lgt1[1]), int(lgt2[1]))
			diff_poffset = math.fabs(int(lgt1[1]) - int(lgt2[1]))
			#print diff_poffset
			
			if x != len(l_poffset) - 1:
				if lowest_poffset in lgt1:
					lowest_poffset_lgt = lgt1
					higher_poffset_lgt = lgt2
				else:
					lowest_poffset_lgt = lgt2
					higher_poffset_lgt = lgt1
			else:
				#print "FIN"
				lowest_poffset_lgt = lgt1
				higher_poffset_lgt = lgt2
			
			print lowest_poffset_lgt
			
			list_metagene_info_lowest_poffset_lgt = []
			
			for g in dict_metagene_profiles[str(lowest_poffset_lgt[0])]:
				list_metagene_info_lowest_poffset_lgt.append(g[1])
				
			print list(reversed(sorted(list_metagene_info_lowest_poffset_lgt)))
			
			list_metagene_info_lowest_poffset_lgt = list(reversed(sorted(list_metagene_info_lowest_poffset_lgt)))
			
			for u in dict_metagene_profiles[str(lowest_poffset_lgt[0])]:
				if u[0] == -lowest_poffset_lgt[1]:
					metagene_info = u[1]
			print metagene_info
			pos_metagene_info = list_metagene_info_lowest_poffset_lgt.index(metagene_info)
			
			new_higher_metagene_info = list_metagene_info_lowest_poffset_lgt[pos_metagene_info + 1]
			print new_higher_metagene_info
			if list_metagene_info_lowest_poffset_lgt.count(new_higher_metagene_info) > 1:
				p_offsets_new_higher_metagene_info = []
			
				print "multiple metagene"
				#diff_offsets_new_higher_metagene_info_previous_poffset = 0
				for e in dict_metagene_profiles[str(lowest_poffset_lgt[0])]:
					if e[1] == new_higher_metagene_info:
						p_offsets_new_higher_metagene_info.append(int(math.fabs(e[0])))
				for t in range(len(p_offsets_new_higher_metagene_info)):
					if t > 0 and t < len(p_offsets_new_higher_metagene_info) - 1:
						current_diff_offsets_new_higher_metagene_info_previous_poffset = p_offsets_new_higher_metagene_info[t]
						next_diff_offsets_new_higher_metagene_info_previous_poffset = p_offsets_new_higher_metagene_info[t + 1]
						
						if current_diff_offsets_new_higher_metagene_info_previous_poffset - higher_poffset_lgt[1] < next_diff_offsets_new_higher_metagene_info_previous_poffset - higher_poffset_lgt[1]:
							new_poffset = current_diff_offsets_new_higher_metagene_info_previous_poffset
						else:
							new_poffset = next_diff_offsets_new_higher_metagene_info_previous_poffset
				print new_poffset
				
				if math.fabs(new_poffset - higher_poffset_lgt[1]) > diff_poffset:
					final_poffset = lowest_poffset[1]
				elif diff_poffset_current_previous < math.fabs(new_poffset - higher_poffset_lgt[1]):		
					final_poffset = lowest_poffset[1]
				else:
					final_poffset = new_poffset
						
			else:
				for e in dict_metagene_profiles[str(lowest_poffset_lgt[0])]:
					if e[1] == new_higher_metagene_info:
						#print new_higher_metagene_info
						new_poffset = int(math.fabs(e[0]))
						#if new_poffset == lowest_poffset_lgt[1]:
						#	new_higher_metagene_info = list_metagene_info_lowest_poffset_lgt[2]
						#	print new_higher_metagene_info
						#	continue
						
						print new_poffset
						if math.fabs(new_poffset - higher_poffset_lgt[1]) > diff_poffset:
							final_poffset = lowest_poffset_lgt[1]
						elif diff_poffset_current_previous != None:
							if diff_poffset_current_previous < math.fabs(new_poffset - higher_poffset_lgt[1]):
								final_poffset = lowest_poffset_lgt[1]
							else:
								final_poffset = new_poffset
						else:
							final_poffset = new_poffset
						print final_poffset	
							#lowest_poffset_lgt = (str(lowest_poffset_lgt[0]), final_poffset)
						dict_final_poffset[str(lowest_poffset_lgt[0])] = final_poffset
					
						if lowest_poffset_lgt == lgt2:
							if x < len(l_poffset) - 1:
								l_poffset[x + 1] = (str(lowest_poffset_lgt[0]), final_poffset)
					
		else:
			print "OK"
			
	print(dict_final_poffset)
	#o.write("length\tp_offset\n")
		
	for k in sorted(dict_final_poffset.keys()):
		o.write(str(k) + "\t" + str(dict_final_poffset[k]) + "\n")
		

	
i.close()
m.close()
o.close()	
			
