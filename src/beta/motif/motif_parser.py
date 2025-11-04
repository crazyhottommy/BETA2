#! /usr/bin/env python
from xml.etree.ElementTree import *
import xml.dom.minidom as minidom
import os, sys
from . import bayesian_motif_comp as bayesian_motif_comp
import numpy
#import mdseqpos.motif as motif

SEP = "|"
def Info(string):
    print(string)

def List2Str(l, con=SEP):
    """use SEP to join list"""
    l = [str(t) for t in l]
    return con.join(l)

def PssmValidator(pssm):
    """validate each PSSM matrix format, no head.
    pssm = [[], [], ... , []]
    """
    #print pssm
    for pos in pssm:
        if len(pos)!=4:
            return False
        for base in pos:
            try:
                float(base)
            except ValueError:
                return False
    return True

def FixPssm(p):
    """FixPssm(p)
    p is a MP object. This function fix the pssm, make each position sum to 1.
    """
    for k in list(p.motifs.keys()):
        for pssm in p.motifs[k]['pssm']:
            for row in pssm:
                rsum = sum([float(t) for t in row])
                for t in range(len(row)):
                    row[t] = round(float(row[t])/rsum,3)
                    if row[t]<0.01:
                        row[t] = 0.01
                imax = row.index(max(row))
                row[imax] -= sum(row)-1

class MotifParser:
    """
    Update 110524:
    for symbol, use :: to seperate components.
    
    MotifParser (MP)
    function: parser motif xml file and return a MP Object "p", you can get the dictionary by "p.motifs"
    If you want to parser a file whose tags are not in pre-decided list, type: p=MotifParser();p.attr_list.append(<sth>);
    
    methods:
    __init__(self, xmlfile = '')	#Setup a MP Object, and parser the xmlfile.
    __add__(self, mp2)				#Add two MP into one, delete duplicates with same motif id. Return a MP object.
    __sub__(self, mp2)				#Subtract motifs of mp2 from the MP object. Return a MP object.
    Parser(self, xmlfile)			#Parser the xmlfile.
    ParserTable(self, tablefile)	#Parser from a "\t" splitted table txt file. Get a MP Object.
    ToXml(self, xmlfile, xsl=False)	#Output the MP Object to xmlfile. if xsl=True, also output the xsl file.
    ToTable(self, tabfile)			#Output the MP Object to table txt file.
    GetAttr(self, attr)             #Get all data of the attribute. Return a list.
    SearchMotif(self, **attrs)      #search a set of specific motifs in the MP Object. Return a MP Object.
    String(self, mid)               #Return a specific motif in a formatted string and readable type.
      
    Details of the dictionary of MP object:
      'id'           list of string    identification in each db.
      'edition'      list of string    edition for JASPAR
      'source'       list of string    from which database
      'sourcefile'   list of string    url link to get the motif in original db
      'status'
      'numseqs'      list of string    length of motif 
      'pmid'         list of string    article support this motif
      'dbd'          list of string    DNA binding domain
      'description'  list of string    gene description
      'species'      list of strings   species the gene belong to
      'entrez'       list of strings   entrezs of the gene
      'symbol'       list of string    symbol of the gene
      'synonym'      list of strings   non-official gene symbols of the gene
      'refseq'       list of strings   refseqs of the gene
      'pssm'         list of 2d-list of float  pssm matrix of the motif
      'comment1'
      'comment2'
      'comment3'
      'comment4'
      'comment5'
    """
    def __init__(self, xmlfile = ''):
        """MP = MotifParser()
        Setup a MP Object, and parser the xmlfile.
        """
        self.motifs = {}
        #self.temp_print = ''
        
        self.keyName = 'id'
        self.attr_list = [self.keyName,'edition']
        self.tag_list = ['source', 'sourcefile', 'status', 'numseqs', 'pmid', 'dbd', \
        'description', 'species', 'cellline', 'entrez', 'symbol', 'synonym', 'refseq', 'comment1', \
        'comment2', 'comment3', 'comment4', 'comment5', 'datasetid', 'zscore', 'seqfactors', \
        'seqdbds', 'nmotifs','logoImg','hits']
        self.special_list = ['pssm'] # if you add a element here, need to edit code below -,-
        self.all_list = self.attr_list + self.tag_list + self.special_list
        
        if xmlfile:
            if xmlfile[-3:] == 'xml':
                self.Parser(xmlfile)
            elif xmlfile[-3:] == 'txt':
                self.ParserTable(xmlfile)
            else:
                Info("Can't parser the file, xml or txt?")
            
    #def Refresh(self):
    #    """Refresh self.all_list, if you ever use sth like self.tag_list.append"""
    #    self.all_list = self.attr_list + self.tag_list + self.special_list

    def AppendTag(self, tag, appendto='tag_list'):
        """AppendTag(self, tag, appendto='tag_list')
        extend tags in MP object.
        """
        taglist = getattr(self, appendto)
        if tag not in taglist:
            taglist.append(tag)
        self.all_list = self.attr_list + self.tag_list + self.special_list
        
    def Create(self, id, **argv):
        """Create(self, id, **argv)
        Create a single motif.
        """
        #self.Refresh()
        self.motifs[id] = {}
        self.motifs[id]['id'] = id
        for i in self.all_list:
            if i in list(argv.keys()):
                self.motifs[id][i] = argv[i][:]
            else:
                self.motifs[id][i] = []
        
    def Parser(self, xmlfile):
        """MP.Parser(xmlfile)
        Parser the xmlfile.
        """
        self.motifs = {}
        xmltree = ElementTree()
        try:
            xmltree.parse(xmlfile)
        except IOError:
            Info("Fail to parser the xml file. Not exist or format error?")
            return None

        for pos in xmltree.findall("motif"):
            #get key and set empty element
            key = pos.get(self.keyName)
            if not key:
                Info('ERROR: No %s found for node.' %self.keyName)
                return None
            if key in list(self.motifs.keys()):
                Info("WARNING: %s has exist in instance."%key)
            self.motifs[key] = {}

            #add attribs and tags for each element.
            for iattr in self.attr_list:
                value = pos.get(iattr)
                if value:
                    self.motifs[key][iattr] = [value]
                else:
                    self.motifs[key][iattr] = []
                
            for itag in self.tag_list:
                self.motifs[key][itag] = [t.text.strip() for t in pos.findall(itag)]
                    
            itag = 'pssm'
            self.motifs[key][itag]=[]
            for ipssm in pos.findall(itag):
                matrix = []
                plist = ipssm.findall('pos')
                plist.sort(key=lambda x:int(x.get('num')))
                for t in plist:
                    base_A = float(t.find('A').text.strip())
                    base_C = float(t.find('C').text.strip())
                    base_G = float(t.find('G').text.strip())
                    base_T = float(t.find('T').text.strip())
                    matrix.append([base_A, base_C, base_G, base_T])
                self.motifs[key][itag].append(matrix)
                        
    def ParserTable(self, tfile):
        """Parser from a "\t" splitted table txt file. 
        The first line should be col name and in lowercase letter.
        The pssm format like this: [[[0.2,0.3,0.3,0.2],[0.1,0.8,0.05,0.05]]]
        """
        self.motifs = {}
        inf = open(tfile)
        line = inf.readline()
        headList = line.rstrip('\n').split('\t')
        headIndex = {}
        for i in range(len(headList)):
            headIndex[headList[i]] = i #headIndex['sourcefile'] = 3
        for each in headList:
            if each not in self.all_list:
                Info("WARNING: column name <%s> is not a node, it will input but can't output. use 'MP_Object.all_list' to get formal format." %each)
                #return 1

        for line in inf:
            linel = line.rstrip('\n').split('\t')
            key = linel[headIndex[self.keyName]] #eg. key = MA00004
            if key in list(self.motifs.keys()):
                Info("WARNING: %s has exist in instance"%key)
            self.motifs[key] = {}

            for iattr in self.attr_list:
                self.motifs[key][iattr] = []
                try:
                    value = linel[headIndex[iattr]]
                    if value:
                        self.motifs[key][iattr] = value.split(SEP)
                except KeyError:
                    pass
            
            for iattr in self.tag_list:
                self.motifs[key][iattr] = []
                try:
                    value = linel[headIndex[iattr]]
                    if value:
                        self.motifs[key][iattr] = value.split(SEP)
                except KeyError:
                    pass

            itag = 'pssm'
            self.motifs[key][itag] = []
            if itag in headList:
                matrix_string = linel[headIndex[itag]]
                if matrix_string:
                    exec('matrix=%s' %matrix_string)
                    for imatrix in matrix:
                        if not PssmValidator(imatrix):
                            Info("ERROR: Matrix format error. id: %s" %key)
                            return False
                    else:
                        self.motifs[key][itag] = matrix
        print("Success parser from table.")
                        
    def GetAttr(self, attr, deldup = False):
        """MP.GetAttr(attr, deldup = False)
        Get all data of the attribute / tag. Return a string.
        """
        if attr not in self.all_list:
            print("Wrong input attr, select attr from: \n:",List2Str(self.all_list, ","))
            return ''
        if attr == 'pssm':
            print("Not support to get pssm. you can use MP.ToTable() ;)\n")
            return ''
        res = []
        for i in list(self.motifs.values()):
            if attr == 'symbol':
                res.extend(i[attr])
            else:
                res.extend(i[attr])
        if deldup:
            res = list(set(res))
        res.sort()
        return List2Str(res).replace("::",SEP)
        
    def SearchMotif(self, **attrs):
        """SearchMotif(self, **attrs)
        search a set of specific motifs in the MP Object. Return a MP Object.
        Any changes on the subset MP may influence the original MP.
        choose arguments from self.all_list
        e.g) MP2 = MP.SearchMotif(species="Homo sapiens",source="JASPAR")
        """
        #print attrs
        for i in list(attrs.keys()):
            if i not in self.all_list:
                Info("Wrong input attr:%s, select attr from:\n: %s" %(i, List2Str(self.attr_list+self.tag_list, ",")))
                return None
        sub_motifs = MotifParser() #self
        sub_motifs.motifs = self.motifs

        for attr in list(attrs.items()):
            temp_dict = {}
            for i in list(sub_motifs.motifs.items()):
                if not attr[1] and not i[1][attr[0]]: #search for empty
                    temp_dict[i[0]] = i[1].copy()
                elif attr[1].upper() in (SEP.join(i[1][attr[0]])).upper().replace('::',SEP).split(SEP):
                    temp_dict[i[0]] = i[1].copy()
            sub_motifs.motifs = temp_dict
        Info("Extract %d records." %len(sub_motifs))
        return sub_motifs

    def String(self, mid):
        """Return a specific motif in a formatted string and readable type.
        e.g) MP.String("M00913")
        """
        if mid in list(self.motifs.keys()):
            dMotif = self.motifs[mid]
        else:
            Info("ID incorrect, can't find Motif ID: %s" %mid)
            return ''
        motif_string = ['\n']
        for itag in self.attr_list + self.tag_list:
            try:
                motif_string.append("%s: %s\n" %(itag, ' '*(10-len(itag)) + List2Str(dMotif[itag]) ))
            except KeyError:
                motif_string.append("%s: None\n" %itag)

        itag = 'pssm'
        for imatrix in dMotif[itag]:
            motif_string.append("PSSM:        A      C      G      T\n")
            for i in range(len(imatrix)):
                motif_string.append("|%6d"%(i+1,) + "  %3.3f  %3.3f  %3.3f  %3.3f\n" %tuple(imatrix[i]))
            motif_string.append("\n")
            
        print(List2Str(motif_string,""))

    def Patch(self, mp2): #marked as not useful
        """patch motif in mp2 to self, simply replace the motif in self, identified by motif id."""
        for item in list(mp2.motifs.items()):
            self.motifs[item[0]] = item[1]
            
    def __add__(self, mp2):
        """MP1.__add__(MP2) <==> MP1+MP2
        Add two MP into one, delete duplicates with same motif id(use MP2 to replace MP1). Return a MP object.
        """
        res_motifs = MotifParser()
        res_motifs.keyName = self.keyName[:]
        res_motifs.attr_list = self.attr_list[:]
        res_motifs.tag_list = self.tag_list[:]
        res_motifs.special_list = self.special_list[:]
        res_motifs.all_list = self.all_list[:]
        
        temp_dict = {}
        for i in list(self.motifs.keys()):
            temp_dict[i] = self.motifs[i].copy()
        for i in list(mp2.motifs.keys()):
            temp_dict[i] = mp2.motifs[i].copy()
        #res_motifs.motifs.update(self.motifs)
        #res_motifs.motifs.update(mp2.motifs)
        res_motifs.motifs = temp_dict
        return res_motifs
    
    def __sub__(self, mp2):
        """MP1.__sub__(MP2) <==> MP1-MP2
        Subtract motifs of mp2 from the MP object (identify by keys). Return a MP object.
        """
        res_motifs = MotifParser()
        res_motifs.keyName = self.keyName[:]
        res_motifs.attr_list = self.attr_list[:]
        res_motifs.tag_list = self.tag_list[:]
        res_motifs.special_list = self.special_list[:]
        res_motifs.all_list = self.all_list[:]

        res_motifs.motifs = self.motifs.copy()
        motif_id_list = list(res_motifs.motifs.keys())
        for i in mp2.motifs:
            if i in motif_id_list:
                del(res_motifs.motifs[i])
        return res_motifs

    def SeqLogo2File(self, folder, rfile = 'draw_seqLogo.r'):
        if not os.path.exists(folder):
            os.mkdir(folder)
        rscript = open(os.path.join(folder, rfile),"w")
        rscript.write('setwd("%s")\n' %folder)
        rscript.write('library("seqLogo")\n')
        for each in list(self.motifs.keys()):
            pssm = self.motifs[each]['pssm']
            if not pssm:
                continue
            if len(pssm)>1:
                print('ERROR: it has more than 1 pssm')
            else:
                pssm = pssm[0]
            t1 = ['c(%s)' %(','.join([str(m) for m in t]),) for t in pssm]
            t2 = 'data<-cbind(%s)\n' %(','.join(t1),)
            rscript.write('png("%s.png", width=660, height=300)\n' %each)
            rscript.write(t2)
            rscript.write('seqLogo(as.matrix(data))\n')
            rscript.write('dev.off()\n\n')

        rscript.close()
        cmd = 'Rscript %s' %os.path.join(folder, rfile)
        os.system(cmd)
        
    def ToXml(self, xmlfile, xsl=False, sortkey=''):
        """MP.ToXML(xmlfile, xsl=False, sortkey='')
        Output the MP Object to xmlfile.
        sortkey = lambda x:x['id'] #something like that.
        """
        doc = minidom.Document()
        motifs = doc.createElement("motifs")
        doc.appendChild(motifs)
        
        t_motifs = list(self.motifs.values())
        if sortkey:
            t_motifs.sort(key=sortkey)
        else:
            #t_motifs.sort(key=lambda x:x[self.keyName]) #sort value as keyName before output
            try:
                t_motifs.sort(key=lambda x:x['symbol'][0].upper()+' '+x['species'][0])
            except:
                t_motifs.sort(key=lambda x:x['symbol'])
            #t_motifs.sort(key=lambda x:x['symbol'][0].upper()+' '+x['species'][0]+' '+'0'*(5-len(x['id'][0][9:]))+x['id'][0][9:]) #sort as symbol, then id, shirley asked
        #t_motifs.sort(key=lambda x:'|'.join(x['dbd'])+' '+x['symbol'][0].upper()) #sort as dbd, then symbol
        for mo in t_motifs:
            # Create the main element
            motif = doc.createElement("motif")
            
            # Create main Attrs
            for iattr in self.attr_list:
                if mo[iattr]:
                    motif.setAttribute(iattr, List2Str(mo[iattr]))
            motifs.appendChild(motif)
    
            #Create elements
            for itag in self.tag_list:
                for ivalue in mo[itag]:
                    element = doc.createElement(itag)
                    motif.appendChild(element)
                    ptext = doc.createTextNode(ivalue)
                    element.appendChild(ptext)
                
            itag = "pssm"
            for matrix in mo[itag]:
                #print matrix
                baseIndex = ['A', 'C', 'G', 'T']
                pssm = doc.createElement(itag)
                motif.appendChild(pssm)
                for i in range(len(matrix)):
                    pos = doc.createElement("pos")
                    pos.setAttribute("num", "%d" %(i+1,))
                    pssm.appendChild(pos)
                    for j in range(4):
                        t = doc.createElement(baseIndex[j])
                        pos.appendChild(t)
                        ptext = doc.createTextNode('%.3f' %matrix[i][j])
                        t.appendChild(ptext)
        xmlString = doc.toprettyxml(indent="  ")
        xmlString = xmlString.split("\n",1)
        
        #output to xmlfile
        xmlfile_ = xmlfile.split(".xml")[0]
        outf = open(xmlfile_+".xml",'w')
        outf.write(xmlString[0]+"\n")
        if xmlString[0].find("xml") != -1:
            outf.write(r'<?xml-stylesheet type="text/xsl" href="%s.xsl"?>'% os.path.split(xmlfile_)[-1] +"\n")
        outf.write(xmlString[1])
        outf.close()
        print("Output xml to file: %s.xml." %xmlfile_)
        
        if xsl:
            outf = open(xmlfile_+".xsl",'w')
            xsl_string_0 = """\
<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="html" encoding="utf-8"/>
<xsl:template match="/">

<html>
<head>
	<title>%s</title>
</head>
<body>
<table cellspacing="1" border="1" bordercolor="#FF9900">
	<THEAD>
		<TR>
			<TD bgcolor="#FFFF99"><B>ID</B></TD>
			<TD bgcolor="#FFFF99"><B>SOURCE</B></TD>
			<TD bgcolor="#FFFF99"><B>STATUS</B></TD>
			<TD bgcolor="#FFFF99"><B>SPECIES</B></TD>
			<TD bgcolor="#FFFF99"><B>CELLLINE</B></TD>
			<TD bgcolor="#FFFF99"><B>SYMBOL</B></TD>
			<TD bgcolor="#FFFF99"><B>ENTREZ</B></TD>
			<TD bgcolor="#FFFF99"><B>REFSEQ</B></TD>
			<TD bgcolor="#FFFF99"><B>DESCRIPTION</B></TD>
			<TD bgcolor="#FFFF99"><B>DBD</B></TD>
			<TD bgcolor="#FFFF99"><B>PMID</B></TD>
			<TD bgcolor="#FFFF99"><B>SEQLOGO</B></TD>
		</TR>
	</THEAD>
	<TBODY>
		<xsl:for-each select="motifs/motif">
		<TR>
			<TD><a href="pwm/{@id}.pwm" target="_blank"><xsl:value-of select="@id" /></a></TD>
			<TD><xsl:if test="not(source)">-</xsl:if>	<xsl:value-of select="source" /></TD>
			<TD><xsl:if test="not(status)">-</xsl:if>	<xsl:value-of select="status" /></TD>
			<TD><xsl:for-each select="species"><xsl:value-of select="." /></xsl:for-each></TD>
			<TD><xsl:for-each select="cellline"><xsl:value-of select="." /></xsl:for-each></TD>
			<TD><xsl:for-each select="symbol"><xsl:value-of select="." /></xsl:for-each></TD>
			<TD><xsl:for-each select="entrez">
			    <a target="_blank"><xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/gene/?term=<xsl:value-of select="normalize-space(.)" />
			    </xsl:attribute><xsl:value-of select="." /></a>
			  </xsl:for-each></TD>
			<TD><xsl:for-each select="refseq"><xsl:value-of select="." /></xsl:for-each></TD>
			<TD><xsl:if test="not(description)">-</xsl:if>	<xsl:value-of select="description" /></TD>
			<TD><xsl:if test="not(dbd)">-</xsl:if>		<xsl:value-of select="dbd" /></TD>
			<TD><xsl:for-each select="pmid">
				<a target="_blank"><xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/pubmed/?term=<xsl:value-of select="normalize-space(.)" />
				</xsl:attribute><xsl:value-of select="." /></a>
			  </xsl:for-each></TD>
			<TD><xsl:if test="not(pssm)">-</xsl:if>		<a href="seqLogo/{@id}.png" target="_blank"><img width="180" height="90" src="seqLogo/{@id}.png"></img></a></TD>
		</TR>
		</xsl:for-each>
	</TBODY>
</table>
</body>
</html>
</xsl:template>
</xsl:stylesheet>"""%xmlfile
            xsl_string_1 = """\
<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="html" encoding="utf-8"/>
<xsl:template match="/">

<html>
<head>
	<title>%s</title>
</head>
<body>
<table cellspacing="1" border="1" bordercolor="#FF9900">
	<THEAD>
		<TR>
			<TD bgcolor="#FFFF99"><B>ID</B></TD>
			<TD bgcolor="#FFFF99"><B>SOURCE</B></TD>
			<TD bgcolor="#FFFF99"><B>STATUS</B></TD>
			<TD bgcolor="#FFFF99"><B>NUMSEQS</B></TD>
			<TD bgcolor="#FFFF99"><B>SPECIES</B></TD>
			<TD bgcolor="#FFFF99"><B>CELLLINE</B></TD>
			<TD bgcolor="#FFFF99"><B>SYMBOL</B></TD>
			<TD bgcolor="#FFFF99"><B>ENTREZ</B></TD>
			<TD bgcolor="#FFFF99"><B>REFSEQ</B></TD>
			<TD bgcolor="#FFFF99"><B>DESCRIPTION</B></TD>
			<TD bgcolor="#FFFF99"><B>PMID</B></TD>
			<TD bgcolor="#FFFF99"><B>RAW_ID</B></TD>
			<TD bgcolor="#FFFF99"><B>SEQLOGO</B></TD>
			<TD bgcolor="#FFFF99"><B>DBD</B></TD>
			<TD bgcolor="#FFFF99"><B>Delete</B></TD>
		</TR>
	</THEAD>
	<TBODY>
		<xsl:for-each select="motifs/motif">
		<TR>
			<TD><xsl:value-of select="@id" /></TD>
			<TD><xsl:if test="not(col2)">-</xsl:if>	<xsl:value-of select="col2" /></TD>
			<TD><xsl:if test="not(status)">-</xsl:if>	<xsl:value-of select="status" /></TD>
			<TD><xsl:if test="not(numseqs)">-</xsl:if>	<xsl:value-of select="numseqs" /></TD>
			<TD><xsl:if test="not(species)">-</xsl:if>	<xsl:for-each select="species"><xsl:value-of select="." /></xsl:for-each></TD>
			<TD><xsl:if test="not(cellline)">-</xsl:if>	<xsl:for-each select="cellline"><xsl:value-of select="." /></xsl:for-each></TD>
			<TD><xsl:if test="not(symbol)">-</xsl:if>	<xsl:value-of select="symbol" /></TD>
			<TD><xsl:if test="not(entrez)">-</xsl:if>	
			  <xsl:for-each select="entrez"><xsl:value-of select="." /></xsl:for-each></TD>
			<TD><xsl:if test="not(refseq)">-</xsl:if>	<xsl:for-each select="refseq"><xsl:value-of select="." /></xsl:for-each></TD>
			<TD><xsl:if test="not(description)">-</xsl:if>	<xsl:value-of select="description" /></TD>
			<TD><xsl:if test="not(pmid)">-</xsl:if><xsl:for-each select="pmid"><xsl:value-of select="." /></xsl:for-each></TD>
			<TD><xsl:if test="not(col1)">-</xsl:if>	<xsl:value-of select="col1" /></TD>
			<TD><xsl:if test="not(pssm)">-</xsl:if>		<a target="_blank"><xsl:attribute name="href">seqLogo/<xsl:value-of select="normalize-space(col1)" />.png</xsl:attribute><img width="180" height="90"><xsl:attribute name="src">seqLogo/<xsl:value-of select="normalize-space(col1)" />.png</xsl:attribute></img></a></TD>
			<TD><xsl:if test="not(dbd)">-</xsl:if>	<xsl:for-each select="dbd"><xsl:value-of select="." /></xsl:for-each></TD>
			<TD><xsl:if test="not(col4)">-</xsl:if>	<xsl:value-of select="col4" /></TD>
		</TR>
		</xsl:for-each>
	</TBODY>
</table>
</body>
</html>
</xsl:template>
</xsl:stylesheet>"""%xmlfile
            outf.write(xsl_string_0)
            outf.close()
            
    def ToTable(self, tabfile, sortkey=''):
        """ToTable(self, tabfile, sortkey='')"""
        pssm_index = self.all_list.index("pssm")
        outf = open(tabfile,'w')
        outf.write(List2Str(self.all_list, "\t")+"\n")
        t_motifs = list(self.motifs.values())
        if sortkey:
            t_motifs.sort(key=sortkey)
        else:
            t_motifs.sort(key=lambda x:x[self.keyName])
        for each in t_motifs:
            motif = [List2Str(each[t]) for t in self.all_list[:pssm_index]]
            if each['pssm']: #output pssm to table
                pssm_slist = []
                for imatrix in each['pssm']:
                    slist = [', '.join(["%.3f"%tt for tt in t]) for t in imatrix]
                    imatrix_string = '[['+'], ['.join(slist)+']]'
                    pssm_slist.append(imatrix_string)
                pssm_string = '['+'], ['.join(pssm_slist)+']'
                motif.append(pssm_string)
            else:
                motif.append('')
            motif.extend([List2Str(each[t]) for t in self.all_list[pssm_index+1:]])
            outf.write(List2Str(motif, "\t")+"\n")
        outf.close()
    
    def Pwm2File(self, folder):
        """Pwm2File(self, folder)
        Output pssm to a folder, filename is the id
        """
        if not os.path.exists(folder):
            os.mkdir(folder)
        for i in list(self.motifs.values()):
            filen = os.path.join(folder, List2Str(i[self.keyName])+'.pwm')
            if len(i['pssm'])!=1:
                Info('Warning: %s may have %d pssm.'%(i[self.keyName][0], len(i['pssm'])))
            if not len(i['pssm']):
                continue
            outf = open(filen, 'w')
            outf.write('A\tC\tG\tT\n')
            for x in i['pssm'][0]:
                x2 = ['%.3d'%t for t in x]
                outf.write(List2Str(x, '\t') + '\n')
            outf.close()
            
    def __len__(self):
        return len(self.motifs)
        
    def _Parser(self, xmlfile):
        """MP.Parser(xmlfile)
        Parser the xmlfile. Old parser version. Only for Convert.
        """
        self.motifs = {}
        xmltree = ElementTree()
        try:
            xmltree.parse(xmlfile)
        except IOError:
            print("Fail to parser the xml file. Not exist or format error?")
            return None

        for pos in xmltree.findall("motif"):
            tag = 'id'
            attrib = pos.attrib
            id = attrib[tag]
            #print id
            self.motifs[id] = {}
            self.motifs[id][tag] = [id]
            
            tag = 'edition'
            try:
                self.motifs[id][tag] = [attrib[tag]]
            except:
                self.motifs[id][tag] = []
    
            for tag in ('source', 'sourcefile', 'status', 'description', 'numseqs', 'pmid', 'dbd', \
            'comment1', 'comment2', 'comment3', 'comment4', 'comment5'):
                if pos.find(tag)!=None:
                    self.motifs[id][tag] = [pos.find(tag).text.strip()]
                else:
                    self.motifs[id][tag] = []
                    
            for tag in ('species', 'entrez', 'synonym', 'refseq', 'symbol'):
                if pos.find(tag+'list')!=None:
                    t = pos.find(tag+'list')
                    self.motifs[id][tag] = [i.text.strip() for i in t.findall(tag)]
                else:
                    self.motifs[id][tag] = []
                    
            for tag in ('pssm',):
                self.motifs[id][tag]=[]
                if pos.find(tag)!=None:
                    plist = pos.find(tag).findall('pos')
                    plist.sort()
                    for t in plist:
                        t = t.findall('*')
                        t.sort()
                        self.motifs[id][tag].append([ float(s.text.strip()) for s in t ])
                    self.motifs[id][tag] = [self.motifs[id][tag]]

    def _ConvertToYing(self):
        """Only For convert."""
        from parser_ying import base_ying
        import numpy
        motiflist = base_ying.MotifList()
        for each in list(self.motifs.values()):
            pm = base_ying.Motif()
            pm.id = each['id'][0]
            pm.status = None
            if each['source']:
                pm.source = each['source'][0]
            if each['sourcefile']:
                pm.sourcefile = each['sourcefile'][0]
            if each['species']:
                pm.species = each['species']
            if each['entrez']:
                pm.entrez = [int(t) for t in each['entrez']]
            if each['symbol']:
                pm.symbol = each['symbol']
            if each['synonym']:
                pm.synonyms = each['synonym']
            if each['description']:
                pm.fullname = each['description'][0]
            if each['dbd']:
                pm.dbd = each['dbd'][0]
            if each['pmid']:
                pm.pmid = each['pmid'][0]
            if each['pssm']:
                pm.pssm = numpy.array(each['pssm'][0], float)
                if pm.pssm.shape[1] != 4:
                    raise ValueError("motif PSSM must have 4 columns")
            pm.numseqs = None
            pm.curators = []
            pm.results = None
            pm.antisense = False # reverse complements need to be taken to make tree consistent
            motiflist.append(pm)
        return motiflist

    def _ConvertToOldMotif(self, motifid):
        d = self.motifs[motifid]
        
        """Returns a new motif with the values specified in the dictionary.
        The dictionary is expected to follow the format specified by the
        MotifParser
        ref: http://cistrome.dfci.harvard.edu/trac/wiki/MotifParser
        """
        fields = ['id', 'status', 'source', 'sourcefile', 'dbd', 'pmid', ]
        lists = ['entrezs', 'refseqs', 'symbols', 'synonyms',]
        #special lists: species, numseqs
        #missing: fullname, curators
        tmp = MotifParser()
        for f in fields:
            #EVERY item in d is a list, so we have to make them into values
            val = ','.join(d[f])
            
            setattr(tmp, f, val)

        for l in lists:
            #we want to set the value, chop off the s for the dictionary
            setattr(tmp, l, d[l[:-1]])


        tmp.pssm = numpy.array(d['pssm'][0])            
        tmp.species = d['species']
        tmp.numseqs = d['numseqs']
        return tmp
    
    def _Similarity(self, motifid1, motifid2, metric='Bayesian'):
        """_Similarity(self, motifid1, motifid2, metric='Bayesian')
        Return a score for the similarity between two motifs.
        offset -- number of basepairs to shift the first motif
        antisense -- whether to take the reverse complement of the first motif
        """
        if len(self.motifs[motifid1]['pssm']) == 1 and len(self.motifs[motifid2]['pssm']) == 1:
            m1 = self._ConvertToOldMotif(motifid1)
            m2 = self._ConvertToOldMotif(motifid2)
            similarity_score, offset, antisense, flag_merge = bayesian_motif_comp.BLiC_score(m1.pssm, m2.pssm)
            antisense = bool(antisense)
            return similarity_score, offset, antisense, flag_merge
        else:
            Info('ERROR: It has no matrix or more than 1 matrix: %s, %s'%(motifid1, motifid2))

    def _ParserTable(self, tfile):
        """Parser from a "\t" splitted table txt file. 
        The first col should be col name and in lowercase letter.
        The pssm format like this: [[0.2,0.3,0.3,0.2],[0.1,0.8,0.05,0.05]]
        """
        self.motifs = {}
        inf = open(tfile)
        line = inf.readline()
        headList = line.rstrip('\n').split('\t')
        headIndex = {}
        for i in range(len(headList)):
            headIndex[headList[i]] = i #headIndex['sourcefile'] = 3
        for each in headList:
            if each not in self.all_list:
                print("WARNING: column name '%s' is not a node, use 'MP_Object.attr_list' to get formal format." %each)
                #return 1
                
        for line in inf:
            linel = line.rstrip('\n').split('\t')
            tag = 'id'
            id = linel[headIndex[tag]]
            self.motifs[id] = {}
            self.motifs[id][tag] = [id]
            
            tag = 'edition'
            self.motifs[id][tag] = []
            if tag in headList:
                if linel[headIndex[tag]]:
                    self.motifs[id][tag] = [linel[headIndex[tag]]] 
                
            for tag in ('source', 'sourcefile', 'status', 'description', 'numseqs', 'pmid', 'dbd', 'symbol', \
            'comment1', 'comment2', 'comment3', 'comment4', 'comment5'):
                self.motifs[id][tag] = []
                if tag in headList:
                    if linel[headIndex[tag]]:
                        self.motifs[id][tag] = [linel[headIndex[tag]]]

            for tag in ('species', 'entrez', 'synonym', 'refseq'):
                self.motifs[id][tag] = []
                if tag in headList:
                    if linel[headIndex[tag]]:
                        self.motifs[id][tag] = linel[headIndex[tag]].split(SEP)
                    
            for tag in ('pssm',):
                pssm = []
                if tag in headList:
                    if linel[headIndex[tag]]:
                        exec('pssm=[%s]' %linel[headIndex[tag]])
                        for ipssm in pssm:
                            if not PssmValidator(ipssm):
                                print("pssm format error. id: %s" %id)
                    self.motifs[id][tag] = pssm

        print("Success parser from table. %s"%tfile)

                
        

