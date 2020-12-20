
/********************************************************************
* Defines a structure for reading strings with command-line option
* style flags -- used for dumping and restoring human-readable
* pedigree data
********************************************************************/

#ifndef FLAGS_H
#define FLAGS_H

#include <unordered_map>
#include <functional>
#include <vector>
#include <string>
#include <list>

#define FLAGS_INPUT_SUCCESS 0
#define FLAGS_INVALID_INPUT 1

typedef std::pair<int, std::function<void(std::vector<std::string>,void*)>> flag_key;

// Structure for reading command-line-flag input and info dumps
struct flag_reader
{
private:
	// Internals
	std::unordered_map<char, std::string> full_names;
	std::unordered_map<std::string, flag_key> exec;
	void* pos = NULL;
public:
	flag_reader() {}
	// Add flag definition
	void add_flag(std::string name, char nick, int narg, std::function<void(std::vector<std::string>,void*)> eff);
	// Read input with flags given
	int read_flags(int narg, char** args);
	int read_flags(std::string args);
	int read_flags(int narg, std::vector<std::string> args);
	// Return if the reader has been initialized
	bool is_new();
	// Possess a flag reader
	void possess(void* p);
	// Get possessor
	void* get_possessor();
};

// Split string on commas
std::vector<std::string> split_opts(std::string s);

#endif
