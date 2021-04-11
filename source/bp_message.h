
/********************************************************************
* Defines the interface provided by a struct representing the
* probability distribution estimates passed as BP 'messages'.
* Overloads arithmetic operators for ease of use.
********************************************************************/

#ifndef BP_MESSAGE_H
#define BP_MESSAGE_H

#include "poisson_pedigree.h"

#include <map>

// Use a canonic form of gene pairs as domain elements of the marginals
// for BP messages
struct bp_domain
{
private:
	std::pair<gene, gene> genes;
public:
	// Constructor
	/// Construct a domain element containing two genes
	bp_domain(gene g1, gene g2);
	// Accessor
	/// Return the ith gene of the domain element
	const gene& operator[](int i) const;
	// Comparators -- for map
	bool operator==(const bp_domain& ot) const;
	bool operator<(const bp_domain& ot) const;
};

// Map that can ignore write requests
template <typename T, typename U>
struct empty_map
{
private:
	std::map<T, U> inner_map;
	bool empty;
public:
	/// Constructors
	empty_map(bool empty);
	empty_map();
	/// Wrap needed operations from map
	std::pair<typename std::map<T, U>::iterator, bool> insert(typename std::pair<T, U> key_val);
	typename std::map<T, U>::iterator find(const T& key);
	typename std::map<T, U>::iterator begin();
	typename std::map<T, U>::iterator end();
	int size() const;
};

// Representation of messages for BP algorithm
struct bp_message
{
private:
	// Private members
	/// Represent non-default values with a map
	empty_map<bp_domain, long double> probabilities;
	/// Represent total probabilities per half-value
	std::map<gene, long double> marginals;
	/// Represent number of unassigned domain elements
	std::map<gene, int> unmapped;
	/// Set default value at initialization time
	long double nullval;
	/// Set domain size (per gene) at initialization time
	int domain_sz;
	// Internal mutation functions
	void update_marginal(const gene& g, const long double& delta);
public:
	// Constructor
	/// Construct a new message with all elements initialized to nullval
	bp_message(long double nullval, int domain_sz);
	// Set & access mapped values
	long double operator[](const bp_domain& value);
	long double get_marginal(const gene& value);
	void inc(const bp_domain& value, const long double delta);
	void set(const bp_domain& value, const long double prob);
	// Arithmetic operations
	/// Add a message to this one
	bp_message& operator+=(bp_message& other);
	/// Multiply this message by another
	bp_message& operator*=(bp_message& other);
	/// Multiply this message by a constant
	bp_message& operator*=(long double scalar);
	/// Divide this message by a constant
	bp_message& operator/=(long double scalar);
	// Normalize this message for future use
	void normalize();
	// Extract maximum value
	bp_domain extract_max();
};

#endif
