#!/usr/bin/env python

# author Firas Said Midnai
# date created  2017-04-11
# date modified 2017-04-17

import numpy as np
import pandas as pd

from collections import OrderedDict

def joinTables(list_of_tables):
    '''
    INPUT
    	list_of_tables : list of pandas.DataFrames
    OUTPUT
    	new_table : pandas.DataFrame of outer concatenation of tables across columns
    EXAMPLE
    	T1_D789 = JoinTables([M[1][df][1] for df in [7,8,9]])
    '''

    new_table = pd.concat(list_of_tables,axis=0,join='inner').fillna(0)
    
    return new_table 

def otuTaxaDict(otu_taxa_map):
	
	'''
	INPUT
		otu_taxa_map : pandas.DataFrame with index of "#OTU ID"s and "taxonomy" column
	OUTPUT
		otu_taxa_dict : pandas.DataFrame where taxonomy has split into distinct levels (p,c,o,f,g,s)	
	NOTES
		* makes sure index for output (i.e. OTU IDs) are integer typed

	EXAMPLE
		otu_taxa_dict = otuTaxaDict(otu_taxa_map)
	'''
	taxa_list = [];
	id_list = [];

	for index,row in otu_taxa_map.iterrows():
	    taxa_dict = {ii.split('__')[0] : ii.split('__')[1] for ii in row['taxonomy'].split('; ')}
	    id_list.append(index)
	    taxa_list.append(taxa_dict)

	otu_taxa_dict = pd.DataFrame(taxa_list,index=id_list)
	otu_taxa_dict.set_index(otu_taxa_dict.index.astype(int),inplace=True)
	otu_taxa_dict.index.name='#OTU ID'
	otu_taxa_dict = otu_taxa_dict.loc[:,['p','c','o','f','g','s']]
	
	return otu_taxa_dict

def pcoa_figure(pcoa,ax,pcs,title,mapping,classes):
    ''' 
    INPUT
        pcoa : dictionary output of read_pcoa_file : {'eigen_vectors' : list, 'eigen_values' : list , 'var_explained' : list]
        ax : figure axis
        pcs : list of string indicating principal components of choice , e.g. ['0','1'] is the first two PCs
        title : string of arbitrary title for plot
        mapping : pandas.DataFrame mappign file
        classes : specifies marker color, shape, and label 
    
    OUTPUT 
        Principal Coordinate Analysis (PCoA) plot, duh ...
    
    '''
    
    class_dict = {};
    
    # grabs coordiantes
    pc_x = pd.Series({k:v[int(pcs[0])] for k,v in pcoa['eigen_vectors'].iteritems()},name='PC'+pcs[0]);
    pc_y = pd.Series({k:v[int(pcs[1])] for k,v in pcoa['eigen_vectors'].iteritems()},name='PC'+pcs[1]);
    
    # tabulate coordinates into pandas.DataFrame
    eigen_vectors = pd.DataFrame([pc_x,pc_y]).T;
    var_explained = pcoa['var_explained'];

    # for each unique group
    for key, value in classes.iteritems():
        # identify data points
        options_dict = value[0]
        #class_dict[key] = mapping[mapping.Column==key].index.values;
        class_dict[key] = mapping[mapping.isin(options_dict).sum(1)==len(options_dict)].index.values
        l_index = [str(ii) for ii in class_dict[key]]
        # plot seperately with unique color, marker, and legend label
        ax.scatter(eigen_vectors.loc[l_index,'PC'+pcs[0]],eigen_vectors.loc[l_index,'PC'+pcs[1]],
                    s=140,color=value[1],edgecolor='black',lw=0.8,label=value[3],zorder=3,marker=value[2])                         

    ax.axhline(0, color='black',alpha=0.5,zorder=1)
    ax.axvline(0, color='black',alpha=0.5,zorder=1)   
    
    x_label_  = 'PC'
    x_label_ += str(int(pcs[0])+1)
    x_label_ += ' ('+ ('%2.2f%%' % float(100*var_explained[int(pcs[0])]))
    x_label_ += ' variation)'

    y_label_  = 'PC '
    y_label_ += str(int(pcs[1])+1)
    y_label_ += ' ('+ ('%2.2f%%' % float(100*var_explained[int(pcs[1])]))
    y_label_ += ' variation)'    
    
    ax.set_xlabel(x_label_,fontsize=20)
    ax.set_ylabel(y_label_,fontsize=20)
    
    ax.set_title(title,fontsize=20,y=1.02,fontweight='bold')

    [tick.label.set_fontsize(14) for tick in ax.xaxis.get_major_ticks()];
    [tick.label.set_fontsize(14) for tick in ax.yaxis.get_major_ticks()];

def phyloSummaryOtuTable(df_parse,otu_taxa_dict,levels=None):
    '''
    DESCRIPTION
    	summary of OTU table based on desired level of phylogenetic binning
    INPUT
    	(1) an OTU table (df_parse)
    	(2) an otu id to phylogeneic mapping (phylum to species)
    	(3) levels of interest, e.g. ['p','g'] for analysis at phylum and genus levels only
        
    OUTPUT
    a dictionary of taxonomy level-specific tables (p,c,o,f,g,s) where 
    each table is clades x the number of unique OTUs and their relative abundances summary statistics
    '''
    df_summary_dict = OrderedDict();

    if type(levels)==list:
    	levels_to_summarize = levels
    elif type(levels)==str:
    	levels_to_summarize = [levels]
    else:
    	levels_to_summarize = ['p','c','o','f','g','s']

    for level in levels_to_summarize:

        if isinstance(df_parse,pd.Series):
                df_summary = pd.DataFrame([],index=otu_taxa_dict.loc[df_parse.keys()].loc[:,level].unique(),
                                 columns=['Number of OTUs','Relative abundance'])
        else:
                df_summary = pd.DataFrame([],index=otu_taxa_dict.loc[df_parse.keys()].loc[:,level].unique(),
                                 columns=['Number of OTUs','mean','median','count','min','max'])

        for ll in otu_taxa_dict.loc[df_parse.keys()].loc[:,level].unique():
            children = otu_taxa_dict[otu_taxa_dict.isin({level:[ll]}).any(1)].index
            children = set(children).intersection(set(df_parse.keys()))
            num_otus = len(children)
            if isinstance(df_parse,pd.Series):
            	df = df_parse.loc[list(children)];
            	df_summary.loc[ll,'Number of OTUs'] = num_otus;
            	df_summary.loc[ll,'Relative abundance'] = df.sum();
            else:
            	df = df_parse.loc[:,list(children)].sum(1)
            	summ_df = [df.mean(0),
                   		   df[df>0].median(0),
                   	       (df>0).sum(0),
                   		   df[df>0].min(0),
                   		   df[df>=0].max(0)];
            	df_summary.loc[ll,:] = [num_otus]+summ_df;
            	

        #print dd,tt,df_summary.shape
        if isinstance(df_parse,pd.Series):
        	df_summary_dict[level] = df_summary.sort_values(['Relative abundance'],ascending=False)
        else:
        	df_summary_dict[level] = df_summary.sort_values(['median'],ascending=False)


    return df_summary_dict

def read_biom(filename):
	'''
	INPUT 
		filename
	OUTPUT
		otu_table: pandas.DataFrame of samples by OTUs
		otu_taxa_map: pandas.DataFrame with index of "#OTU ID"s and "taxonomy" column
	EXAMPLE
		otu_table,otu_taxa_map = read_biom('./tables/otus_table.filtered.from_biom.txt')
	'''
	biom = pd.read_csv(filename,sep='\t',header=0,index_col=0,skiprows=1);
	otu_taxa_map = pd.DataFrame(biom.iloc[:,-1]);
	otu_table = biom.iloc[:,:-1];

	return otu_table,otu_taxa_map

def read_pcoa_file(filename):
	'''
	INPUT EXAMPLE
		"./pc/pcoa_binary_jaccard_otus_table.filtered.txt"
	OUTPUT
		a dictionary with three keys
			eigen_vectors is a dictionary where keys are samples and values are lists corresponding to the coordinates of a sample on a principal comopnents (ordered first to last)
			eigen_values is a list of the eigen values of the principal components (from first to last)
			var_explained is a list of the cumulative variances explained by principal components (from first through last)
	'''

	fid = open(filename,'r')

	_eigvec_  = 0;
	_eigval_  = 0;
	_explain_ = 0;

	eigen_vectors  = {};
	eigvals        = [];
	var_explained  = [];
	pcoa_dict      = {};

	for line in fid:

	    if _eigvec_==1:

	        if line =='\n':
	            _eigvec_=0;
	            continue;

	        line       = line.strip()
	        line_key   = line.split('\t')[0]
	        line_value = [float(ii) for ii in line.split('\t')[1:]];

	        eigen_vectors.update({line_key:line_value});

	    elif _eigval_==1:

	        eigen_values  = [float(ii) for ii in line.strip().split('\t')]; 
	        _eigval_      = 0;

	    elif _explain_==1:

	        var_explained  = [float(ii) for ii in line.strip().split('\t')]; 
	        _explain_      = 0;

	    #endif

	    if   line.split('\t')[0]=="Site":

	         _eigvec_ = 1;        

	    elif line.split('\t')[0]=='Eigvals':

	         _eigval_ = 1;

	    elif line.split('\t')[0]=='Proportion explained':

	         _explain_ = 1;

	    #endif

	pcoa_dict['eigen_values']  = eigen_values;
	pcoa_dict['eigen_vectors'] = eigen_vectors;
	pcoa_dict['var_explained'] = var_explained;

	return pcoa_dict

def subsetTableByMetadata(otu_df,mapping_df,options_dict):
	
	'''
	INPUT 
		otu_df : pandas.DataFrame of samples x OTUs
		mapping_df : pandas.DataFrame with mapping of attributes in "options_dict"
		options_dict : dictionary where keys are columns in mappping_df and values are lists of column values of interest
	OUTPUT
		microcosms_otu_df : pandas.DataFrame which is the whole or subset of otu_df based on criteria determined in "options_dict"

	EXAMPLE
		subsetTableByMetadata(biom,mapping,{'Time/Comment':['0'],'Column':[1]})
	'''

	microcosms = mapping_df[mapping_df.isin(options_dict).sum(1)==len(options_dict)]
	microcosms_otu_df = otu_df.loc[[str(ii) for ii in microcosms.index],:]

	# below removes OTUs with zero counts
	todrop = microcosms_otu_df.columns[np.where(microcosms_otu_df.sum()==0)[0]]
	microcosms_otu_df = microcosms_otu_df.drop(todrop,axis=1)

	return microcosms_otu_df

def subsetTableBySampleIDs(otu_df,sample_ids):

	'''
	INPUT
		otu_df : pandas.DataFrame of samples x OTUs
		sample_ids : list of sample IDs
	OTUPUT
		otu_df : pandas.DataFrame which is the whole or subset of otu_df based on sample ids selected

	EXAMPLE
		subsetTableBySampleIDs(T0_D1,['80','111'])
	'''

	if len(set(otu_df.index).intersection(set(sample_ids)))!=len(sample_ids):
		return 'ERROR : Some samples are missing from the table.'
	else: 
		otu_df = otu_df.loc[sample_ids,:]

		# below removes OTUs with zero counts
		todrop = otu_df.columns[np.where(otu_df.sum()==0)[0]]
		otu_df = otu_df.drop(todrop,axis=1)
		
		return otu_df.loc[sample_ids,:]

def summarizeTable(df,otu_taxa_dict):
	'''
	INPUT
		df : pandas.DataFrame of samples by OTUs
		otu_taxa_dict : pandas.DataFrame where taxonomy has split into distinct levels (p,c,o,f,g,s)
	OUTPUT
	EXAMPLE
	'''
	summ_df = pd.DataFrame([df.mean(0),
							df[df>0].median(0),
							df[df>0].min(0),
							df[df>=0].max(0),
							(df>0).sum(0)],
							index=['mean','median','min','max','count']).T;
	summ_df = summ_df.join(otu_taxa_dict);
	summ_df = summ_df.loc[:,['mean','median','count','min','max','p','c','o','f','g','s']];
	return summ_df

def tss_norm(otu_table,axis=1):

    '''
    INPUT
        otu_table: pandas.DataFramE
        axis: either 0 or 1. 0 assumes that that otu_table is formatted by samples x otus and therefore sums across each row. 
    OUTPUT
        otu_table: pandas.DataFrame where rows (or columns) sum to one. 
    EXAMPLE
        normalized_otu_table = tss_norm(otu_table,1)
    NOTES
        tss stands for Total Sum Scaling
    '''

    # sum across each row
    if axis==1:
        new_table = otu_table.T / otu_table.sum(1);
    # sum across each column
    elif axis==0:
        new_table = otu_table / otu_table.sum(0);
    else: 
        return 'ERROR: axis can only be 0 or 1.'
    return new_table.T
