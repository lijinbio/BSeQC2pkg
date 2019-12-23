#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <iterator>
#include <thread>
#include <cstdio>
#include "sam.h"

using namespace boost::program_options;
using namespace std;

class Opts {
	public:
		string infile;
		string outfile;
		string reference;
		int length;
		int numthreads;
	public:
		Opts():
			infile("")
			, outfile("")
			, reference("")
			, length(150)
			, numthreads(1) { }
	public:
		void out() {
			cout << "infile: " << infile << endl;
			cout << "outfile: " << outfile << endl;
			cout << "reference: " << reference << endl;
			cout << "length: " << length << endl;
			cout << "numthreads: " << numthreads << endl;
		}
} opts;

int parse_options(int ac, const char ** av) {
	try
	{
		options_description desc{"Allowed options"};
		desc.add_options()
			("help,h", "Produce help message.")
			("infile,i", value<string>()->default_value(""), "Input BAM file. It should be indexed.")
			("outfile,o", value<string>()->default_value(""), "Output statistics.")
			("reference,r", value<string>()->default_value(""), "Reference FASTA file.")
			("length,l", value<int>()->default_value(150), "Read length. Default: 150.")
			("numthreads,t", value<int>()->default_value(1), "Number of threads. Default: 1.")
			;

		variables_map vm;
		store(parse_command_line(ac, av, desc), vm);
		notify(vm);

		if (vm.count("help")) {
			cout << desc << endl;
			cout << "Examples:" <<endl;
			cout << "  " << av[0] << " -i in.bam -o readcount.txt -t 8 -l 160 -r hg38.fa" << endl;
			cout << endl;
			cout << "Date: 2019/12/17" << endl;
			cout << "Authors: Jin Li <lijin.abc@gmail.com>" << endl;
			exit(1);
		}

		for(map<string, variable_value>::iterator it=vm.begin(); it!=vm.end(); ++it) {
			string k=it->first;
			if(k=="infile"){
				opts.infile=vm[k].as<string>();
			} else if(k=="length"){
				opts.length=vm[k].as<int>();
			} else if(k=="numthreads"){
				opts.numthreads=vm[k].as<int>();
			} else if(k=="reference"){
				opts.reference=vm[k].as<string>();
			} else if(k=="outfile"){
				opts.outfile=vm[k].as<string>();
			} else {
				cerr << "Error: invalid option " << k << endl;
				exit(1);
			}
		}
		if (opts.infile.empty()) {
			cerr << "Error: -i|--infile must be specified." << endl;
			cout << desc << endl;
			exit(1);
		}
		if (opts.outfile.empty()) {
			cerr << "Error: -o|--outfile must be specified." << endl;
			cout << desc << endl;
			exit(1);
		}
		if (opts.reference.empty()) {
			cerr << "Error: -r|--reference must be specified." << endl;
			cout << desc << endl;
			exit(1);
		}
		opts.out();
	} catch (const error &ex) {
		cerr << ex.what() << endl;
	}
	return 0;
}

void refbychr(string infile, string chr, string & ref)
{
	ifstream fin(infile);
	string line;
	while ((fin.good() && !fin.eof()) && getline(fin, line)) {
		if (line.empty()) continue;
		if (line==">"+chr) {
			while ((fin.good() && !fin.eof()) && getline(fin, line)) {
				if (line.empty()) continue;
				if (line[0]=='>') break;
				ref+=line;
			}
		}
	}
	fin.close();
}

// Global map for each chr, remember to clear for each chr
map< int, map< string, vector< vector< int > > > > chrtagreadcounts; // chrid->tag->[[c1,...,clength],[t1,...,tlength]]
static int addtag(const bam1_t *b, void *data) {
	if(b->core.n_cigar > 3) return 1;

	string &ref = *(static_cast< string* >(data));
	size_t qlen=b->core.l_qseq;
	if (opts.length<qlen) {
		cerr << "Error: the input length is less than the read length " << qlen << endl;
		exit(1);
	}

	char * s = (char*) bam1_seq(b);
	string seq(qlen, '\0');
	for (int i = 0; i<qlen; ++i) {
		seq[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
	}

	string newSeq;
	int refLen=0;
	int spos=0;
	int okToProceed=1;
	uint32_t* cigar=bam1_cigar(b);
	for( int k=0; k< b->core.n_cigar; ++k)
	{
		int cop =cigar[k] & BAM_CIGAR_MASK; // operation
		int cl = cigar[k] >> BAM_CIGAR_SHIFT; // length
		switch(cop)
		{
			case BAM_CMATCH:
				newSeq  += seq.substr(spos, cl );
				spos += cl;
				refLen += cl;
				break;
			case BAM_CINS:
				spos += cl;
				break;
			case BAM_CDEL:
				for(int itemp = 0; itemp < cl; itemp++) { newSeq += 'N'; }
				refLen += cl;
				break;
			case BAM_CREF_SKIP:
				okToProceed = 0;
				break;
			case BAM_CSOFT_CLIP:
				okToProceed = 0;
				break;
			case BAM_CHARD_CLIP:
				okToProceed = 0;
				break;
			case BAM_CPAD:
				okToProceed = 0;
				break;
			default:
				okToProceed = 0;
				break;
		}
	}

	if(! okToProceed) return 1;
	seq=newSeq;
	string tag=string((char *)bam_aux2Z(bam_aux_get(b, "ZS")));
	string refseq;
	bool forwardstrand=(tag[0]=='+');
	bool forwardsave=(tag=="++" || tag=="--");
	if (forwardstrand) {
		refseq=ref.substr(b->core.pos, refLen+1);
	} else {
		if (b->core.pos==0) {
			refseq="N"+ref.substr(0, refLen);
		} else {
			refseq=ref.substr(b->core.pos-1, refLen+1);
		}
	}
	boost::to_upper(refseq);

	map< string, vector< vector< int > > > &tagreadcounts=chrtagreadcounts[b->core.tid]; // tag->[[c1,...,clength],[t1,...,tlength]]
	map< string, vector< vector< int > > > :: iterator it=tagreadcounts.find(tag);
	if (tagreadcounts.end()!=it) {
		if (forwardstrand) {
			if (forwardsave) {
				for (int i=0; i<qlen; ++i) {
					if (refseq[i]=='C') {
						if (refseq[i+1]=='G') {
							it->second[0][i]+=(seq[i]=='C');
							it->second[1][i]+=(seq[i]=='T');
						} else {
							it->second[2][i]+=(seq[i]=='C');
							it->second[3][i]+=(seq[i]=='T');
						}
					}
				}
			} else {
				for (int i=0; i<qlen; ++i) {
					if (refseq[i]=='C') {
						if (refseq[i+1]=='G') {
							it->second[0][qlen-1-i]+=(seq[i]=='C');
							it->second[1][qlen-1-i]+=(seq[i]=='T');
						} else {
							it->second[2][qlen-1-i]+=(seq[i]=='C');
							it->second[3][qlen-1-i]+=(seq[i]=='T');
						}
					}
				}
			}
		} else {
			if (forwardsave) {
				for (int i=0; i<qlen; ++i) {
					if (refseq[i+1]=='G') {
						if (refseq[i]=='C') {
							it->second[0][i]+=(seq[i]=='G');
							it->second[1][i]+=(seq[i]=='A');
						} else {
							it->second[2][i]+=(seq[i]=='G');
							it->second[3][i]+=(seq[i]=='A');
						}
					}
				}
			} else {
				for (int i=0; i<qlen; ++i) {
					if (refseq[i+1]=='G') {
						if (refseq[i]=='C') {
							it->second[0][qlen-1-i]+=(seq[i]=='G');
							it->second[1][qlen-1-i]+=(seq[i]=='A');
						} else {
							it->second[2][qlen-1-i]+=(seq[i]=='G');
							it->second[3][qlen-1-i]+=(seq[i]=='A');
						}
					}
				}
			}
		}
	} else {
		vector< vector< int > > readcounts(4, vector< int > (opts.length, 0));
		if (forwardstrand) {
			if (forwardsave) {
				for (int i=0; i<qlen; ++i) {
					if (refseq[i]=='C') {
						if (refseq[i+1]=='G') {
							readcounts[0][i]+=(seq[i]=='C'); // CG
							readcounts[1][i]+=(seq[i]=='T');
						} else {
							readcounts[2][i]+=(seq[i]=='C'); // CH
							readcounts[3][i]+=(seq[i]=='T');
						}
					}
				}
			} else {
				for (int i=0; i<qlen; ++i) {
					if (refseq[i]=='C') {
						if (refseq[i+1]=='G') {
							readcounts[0][qlen-1-i]+=(seq[i]=='C'); // CG
							readcounts[1][qlen-1-i]+=(seq[i]=='T');
						} else {
							readcounts[2][qlen-1-i]+=(seq[i]=='C'); // CH
							readcounts[3][qlen-1-i]+=(seq[i]=='T');
						}
					}
				}
			}
		} else {
			if (forwardsave) {
				for (int i=0; i<qlen; ++i) {
					if (refseq[i+1]=='G') {
						if (refseq[i]=='C') {
							readcounts[0][i]+=(seq[i]=='G'); // CG
							readcounts[1][i]+=(seq[i]=='A');
						} else {
							readcounts[2][i]+=(seq[i]=='G'); // CH
							readcounts[3][i]+=(seq[i]=='A');
						}
					}
				}
			} else {
				for (int i=0; i<qlen; ++i) {
					if (refseq[i+1]=='G') {
						if (refseq[i]=='C') {
							readcounts[0][qlen-1-i]+=(seq[i]=='G'); // CG
							readcounts[1][qlen-1-i]+=(seq[i]=='A');
						} else {
							readcounts[2][qlen-1-i]+=(seq[i]=='G'); // CH
							readcounts[3][qlen-1-i]+=(seq[i]=='A');
						}
					}
				}
			}
		}
		tagreadcounts[tag]=readcounts;
	}
	return 0;
}

void bseqcreadcountschrbatch(string bamfile, vector< string > chrs) {
	samfile_t *in=0;
	if ((in=samopen(bamfile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << bamfile << endl;
		return;
	}

	map< string, int > chr2tid;
	for (int i=0; i<in->header->n_targets; i++) {
		chr2tid[in->header->target_name[i]]=i;
	}

	bam_index_t *idx=0;
	idx = bam_index_load(bamfile.c_str());
	if (idx==0) {
		cerr << "Error: not found index file of " << bamfile << endl;
		return;
	}
	for (string &chr : chrs) {
		cout << "Start chromosome " << chr << endl;
		string ref;
		refbychr(opts.reference, chr, ref);
		int tid, beg, end, result;
		bam_parse_region(in->header, chr.c_str(), &tid, &beg, &end);
		if (tid<0) { 
			cerr << "Error: unknown reference name " << chr << endl;
			return;
		}
		result=bam_fetch(in->x.bam, idx, tid, beg, end, static_cast<void*>(&ref), addtag);
		if (result<0) {
			cerr << "Error: failed to retrieve region " << chr << endl;
			return;
		}
		cout << "End chromosome " << chr << endl;
	}
	samclose(in);
	bam_index_destroy(idx);
}

int bseqcreadcounts(string bamfile)
{
	samfile_t *in=0;
	if ((in=samopen(bamfile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << bamfile << endl;
		return 1;
	}
	vector< string> chroms;
	for (int i=0; i<in->header->n_targets; i++) {
		chroms.push_back(in->header->target_name[i]);
	}
	samclose(in);

	vector< vector< string > > chrbatch;
	for (int i=0; i<chroms.size(); i++) {
		if (i>=opts.numthreads) {
			chrbatch[i%opts.numthreads].push_back(chroms[i]);
		} else {
			vector< string > chrs {chroms[i]};
			chrbatch.push_back(chrs);
		}
	}

	vector<thread> threads;
	for (vector< string > &chrs : chrbatch) {
		threads.push_back(thread(bseqcreadcountschrbatch, bamfile, chrs));
	}
	for (auto& th : threads) {
		th.join();
	}

	map< string, vector< vector< int > > > tagsresult;
	for (map< int, map< string, vector< vector< int > > > > :: iterator it=chrtagreadcounts.begin(); chrtagreadcounts.end()!=it; ++it) {
		map< string, vector< vector< int > > > &tagreadcounts=it->second;
		for (map< string, vector< vector< int > > > :: iterator ittag=tagreadcounts.begin(); tagreadcounts.end()!=ittag; ++ittag) {
			string tag=ittag->first;
			vector< vector< int > > &rctag=ittag->second;

			map< string, vector< vector< int > > > :: iterator ithit=tagsresult.find(tag);
			if (tagsresult.end()!=ithit) {
				vector< vector< int > > &readcounts=ithit->second;
				for (int i=0; i<readcounts.size(); i++) {
					for (int j=0; j<readcounts[i].size(); j++) {
						readcounts[i][j]+=rctag[i][j];
					}
				}
			} else {
				tagsresult[tag]=rctag;
			}
		}
	}
	ofstream fout(opts.outfile);
	fout << "position" << '\t' << "numC_CG" << '\t' << "numT_CG" << '\t' << "numC_CH" << '\t' << "numT_CH" << '\t' << "tag" << endl;
	for (map< string, vector< vector< int > > > :: iterator it=tagsresult.begin(); it!=tagsresult.end(); ++it) {
		string tag=it->first;
		vector< int > &ccg=it->second[0];
		vector< int > &tcg=it->second[1];
		vector< int > &cch=it->second[2];
		vector< int > &tch=it->second[3];
		for (int i=0; i<ccg.size(); i++) {
			fout << i+1 << '\t' << ccg[i] << '\t' << tcg[i] << '\t' << cch[i] << '\t' << tch[i] << '\t' << tag << endl;
		}
	}
	return 0;
}

int main(int argc, const char ** argv)
{
	parse_options(argc, argv);
	bseqcreadcounts(opts.infile);
	return 0;
}
