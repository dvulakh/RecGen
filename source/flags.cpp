
/********************************************************************
* Implements a structure for reading strings with command-line
* option style flags -- used for dumping and restoring human-
* readable pedigree data
********************************************************************/

#include "flags.h"
#include <sstream>

// Add new flag to nickname map and effects map
void flag_reader::add_flag(std::string name, char nick, int narg, std::function<void(std::vector<std::string>,void*)> eff)
{
    exec.insert({ name, { narg, eff } });
    if (nick)
        full_names.insert({ nick, name });
}

// Read input string
int flag_reader::read_flags(int narg, std::vector<std::string> args)
{
    std::list<std::string> flag_queue;
    int sn = 0;
    while (sn < narg) {
        std::string s = args[sn];
        /// New flag(s)
        if (s[0] == '-') {
            /// No flags
            if (s.length() < 2)
                return FLAGS_INVALID_INPUT;
            /// Long flag
            else if (s[1] == '-')
                flag_queue.push_back(s.substr(2, s.length()));
            /// Small flags
            else
                for (int i = 1; i < s.length(); i++) {
                    auto it = full_names.find(s[i]);
                    if (it == full_names.end())
                        return FLAGS_INVALID_INPUT;
                    flag_queue.push_back(it->second);
                }
            sn++;
        }
        /// New argument(s)
        else {
            /// No flags waiting -- pass instead of raising errors
            if (flag_queue.empty()) {
                sn++;
                continue;
            }
            /// Execute function
            std::string f = flag_queue.front(); flag_queue.pop_front();
            auto it = exec.find(f);
            if (it == exec.end())
                return FLAGS_INVALID_INPUT;
            flag_key fk = it->second;
            std::vector<std::string> ar;
            int n_tok = fk.first;
            if (n_tok < 0)
                try { n_tok = std::stoi(args[sn++]); }
                catch (std::exception e) { return FLAGS_INVALID_INPUT; }
            for (int i = 0; i < n_tok; i++)
                if (sn + i >= narg)
                    return FLAGS_INVALID_INPUT;
                else
                    ar.push_back(args[sn + i]);
            try { fk.second(ar, pos); }
            catch (std::exception e) { return FLAGS_INVALID_INPUT; }
            sn += n_tok;
        }
    }
    /// Leftover flags must take no arguments
    for (std::string s : flag_queue) {
        auto it = exec.find(s);
        if (it == exec.end() || it->second.first)
            return FLAGS_INVALID_INPUT;
        try { it->second.second(std::vector<std::string>(), pos); }
        catch (std::exception e) { return FLAGS_INVALID_INPUT; }
    }
    return FLAGS_INPUT_SUCCESS;
}

// Make a vector and send to the reader
int flag_reader::read_flags(int narg, char** args)
{
    std::vector<std::string> argv;
    for (int i = 1; i < narg; i++)
        argv.push_back(std::string(args[i]));
    return this->read_flags(narg - 1, argv);
}

// Make a vector and send to the reader
int flag_reader::read_flags(std::string args)
{
    std::vector<std::string> argv;
    std::istringstream iss(args);
    for(std::string s; iss >> s;)
        argv.push_back(s);
    return this->read_flags(argv.size(), argv);
}

// Return if reader has been initialized
bool flag_reader::is_new() { return this->exec.empty(); }

// Possess a flag reader
void flag_reader::possess(void* p) { this->pos = p; }

// Get possessor
void* flag_reader::get_possessor() { return this->pos; }

// Split string on commas
std::vector<std::string> split_opts(std::string s)
{
    std::stringstream sin(s);
    std::vector<std::string> opts;
    while (sin.good()) {
        std::string o;
        getline(sin, o, ',');
        opts.push_back(o);
    }
    return opts;
}
