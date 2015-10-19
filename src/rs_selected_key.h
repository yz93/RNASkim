#ifndef RS_SELECTED_KEY_H
#define RS_SELECTED_KEY_H
#include <string>
#include <vector>
#include <map>

using std::endl;
using std::fstream;
using std::ios;
using std::map;
using std::string;
using std::vector;

namespace rs {
	class TranscriptInfo{
	public:
		TranscriptInfo(int tidx_param) : tidx(tidx_param) {}
		int tidx;  // transcript index
		vector<int> positions;  // the owner k-mer's positions in tidx
	};


	// Key represents a k-mer class
	class Key{
	public:
		Key(string key_param, int count_param) : key(key_param), count(count_param) {}
		string key;  // k-mer sequence
		vector<TranscriptInfo> transcript_infos;  // its associated transcripts information
		int count;  // count of this k-mer (in a gene)
		int transcript_infos_size() const { return transcript_infos.size(); }
	};

	class SelectedKey{
	public:
		string gid;  // Gene name
		vector<Key> keys;
		vector<string> tids;
		vector<int> lengths;  // records the length of corresponding (at same array index) transcript in tids.
		int tids_size() const { return tids.size(); }
		int keys_size() const { return keys.size(); }
	};

	map<string, SelectedKey> process_selected_keys(string cluster_file, string sigmer_count_file, string sigmer_details);
}// namespace rs


#endif