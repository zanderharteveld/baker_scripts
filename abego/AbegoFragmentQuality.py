#!/usr/bin/env  python
import sys,os,glob,re
from Blueprint import Blueprint
import argparse
import string
from operator import add

def pp_class(ppo): ## E G A B O
    ppo = ( In_range( ppo[0]), In_range(ppo[1]), In_range(ppo[2]))
    assert -180<=ppo[0]<=180 and -180<=ppo[1]<=180 and -180<=ppo[1]<=180

    if abs( ppo[2] ) < 90: return 'O'
    elif ppo[0]>=0:
        if -100< ppo[1] <= 100:return 'G'
        else: return 'E'
    elif -125 < ppo[1] <= 50: return 'A'
    else: return 'B'

def In_range(angle):
    while angle>180: angle = angle-360
    while angle<=-180:angle = angle+360
    return angle

def MakeBlueprint(pdb):
        newbluefile = pdb + '.bp'
        cmd = "/work/embqtc/rosetta/rosetta_source/bin/make_blueprint.static.linuxiccrelease -s %s  -database /work/embqtc/rosetta/rosetta_database/ -blue %s -mute all" %( pdb, newbluefile )
	os.system( cmd )

def MakeAbegoList(newbluefile):
        blue = Blueprint( newbluefile )
	blue.reindex_blueprint(start=1)
        abego_list=[]
        for i,res in enumerate( blue.bp_data ):
                abego_list.append( res[2][1] )
        return abego_list


def AbegoFragQuality( pdb, fragment_file, outfilename,**kwargs):
	bluefile = kwargs.get('bluefile')
	if bluefile == None:
		MakeBlueprint(pdb)
		bluefile = pdb+'.bp'

        abego_list = MakeAbegoList(bluefile)

        # frag_ss: list of 3 objects (e,h, l). Each member corresponds to a list of probabilities for all positions to be a given SS (in order: E, H, L)
        # abgeo: the same for abegos in the following order: A, B, G, E, O

        frag_ss,abgeo_frags = Read_fragments(fragment_file)
        #--------------
        abego_types = 'ABGEO'
        probabilities = []
        outfile = open(outfilename,'w')
        for i,abego in enumerate(abego_list):
                index = abego_types.find(abego)
                prob = abgeo_frags[index][i]
                outfile.write('%i %s %s\n' %( i+1, abego, prob))
        probabilities.append(prob)
        outfile.close()

def Read_fragments(fragment_file):
# 1st five lines of frag file:
# position: 1 neighbors: 200
#
# 1di1 A 229 S L -147.015 136.320 177.891
# 1di1 A 230 A H -61.103 -28.146 179.117
# 1di1 A 231 V H -66.832 -46.265 179.757
    FRAG_SIZE = 9
    base = string.split(fragment_file,'/')[-1]
    size = FRAG_SIZE #int(base[7:9])
    sys.stderr.write('Reading fragment file: %s size= %d\n'\
                     %(fragment_file,size))
    if fragment_file[-3:] == '.gz':
        data = popen('gzip -dc '+fragment_file)
    else:
        data = open(fragment_file,'r')
    line = data.readline()
    prev = (-1,-1)
    ss_count = {}
    ppo_count = {}
    while line:
        l = string.split(line)
        line = data.readline()
        assert len(l) == 4 and l[0] == 'position:'

        window = int(l[1])-1 ## numbering starts at 0
        nbrs = int(l[3])
        for i in range(size):
            pos = window + i
            if not ss_count.has_key( pos ):
                ss_count[ pos ] = {'H':0,'E':0,'L':0}
                ppo_count[ pos ] = {'A':0,'B':0,'G':0,'E':0,'O':0}

        for n in range(nbrs):
            for i in range(size):
                line = data.readline()
                l = string.split(line)
                pos = window + i
                ss = l[4]
                ppo = pp_class( map(float, l[5:8] ) )
                ss_count[pos][ss] += 1
                ppo_count[pos][ppo] += 1
            line = data.readline() ## extra blank line

        line = data.readline()
    data.close()

    L = len(ss_count.keys())

    l = []
    e = []
    h = []
    abgeo = ( [], [], [], [], [] )
    ABGEO = 'ABGEO'
    for i in range(L):
        total = reduce(add,ss_count[i].values())
        if total>0:
            e.append( float(ss_count[i]['E'])/total)
            h.append( float(ss_count[i]['H'])/total)
            l.append( float(ss_count[i]['L'])/total)
            for j in range(len(ABGEO)):
                count = ppo_count[i][ABGEO[j]]
                frac = float( count )/total
                if not (frac>= 0 and frac<=1):
                    print i,count,total,frac
                abgeo[j].append( frac )

    frag_ss = (e,h,l)

    return frag_ss,abgeo

def ProbabilityAbegoFragments(abego_frag_qual, start,end):
        lines = open(abego_frag_qual).readlines()
        loop_lines = lines[start:end+1]
        probs=[] ; total_p = 1.0
        for line in loop_lines:
                p = float(line.split()[2])
                total_p *= p
                probs.append(p)

        return total_p


parser = argparse.ArgumentParser()
parser.add_argument('-pdb', type=str)
parser.add_argument('-frag9', type=str)
parser.add_argument('-output', type=str)

args = parser.parse_args()

AbegoFragQuality( args.pdb, args.frag9, args.output)

