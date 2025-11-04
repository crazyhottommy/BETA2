def get_length(inf):
    
    chrom_locate = {}
    count = 0
    linum = 0
    
    for line in inf:
        linelen = len(line)
        count += linelen
        linum += 1
        if line.startswith('>'):
            chrom_locate[line] = [linum,count]
    
    b = sorted(iter(chrom_locate.items()), key=lambda d:d[0], reverse = False)#the sorted chrom_locate, is a list
    for i in range(len(b)-1):
        totalline = b[i+1][1][0] - b[i][1][0] - 1#how many lines of this chromsome
        totalbase = b[i+1][1][1] - b[i][1][1] - len(b[i][0])#how many bases in this chromsome
        width = totalbase/totalline
        b[i][1].append(width)
    last_width = (count - b[len(b)-1][1][1])/(linum - b[len(b)-1][1][0])
    b[len(b)-1][1].append(last_width)
    return b
    #b = [('>chr1\n',[line,base,width]),(),()]
def seeker(inf,chrom_locate,chrom,start,end):
    region = int(end) - int(start)
    chroms = '>' + chrom + '\n'
    sequence = ''
    for i in chrom_locate:
        if i[0] == chroms:
            fs = i[1][1]#fs:first start, the first site of the sequence line
            width = i[1][-1]
            n_before_region = int(start)/width
            st = int(fs) + int(start) + n_before_region#start
            inf.seek(st)#start from st
            sequence = inf.read(region*2).replace('\n','')[0:region]
        else:
            continue
    return sequence

def runfastafrombed(genome, regionf, output):
    
    inf = open(genome)
    chrom_locate = get_length(inf)
    bedf = open(regionf)
    outf = open(output,'w')
    for line in bedf:
        line = line.strip()
        line = line.split('\t')
        chrom = line[0]
        start = line[1]
        end = line[2]
        sequence = seeker(inf,chrom_locate,chrom,start,end)
        if sequence != '':
        	outf.write('>%s:%s-%s\n'%(chrom,start,end))
        	outf.write(sequence + '\n')
        else:
        	pass
    outf.close()
    
