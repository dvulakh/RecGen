
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
	const gene& operator[](int i);
	// Comparators -- for map
	bool operator==(const bp_domain& ot) const;
	bool operator<(const bp_domain& ot) const;
};

// Representation of messages for BP algorithm
struct bp_message
{
private:
	/// Represent non-default values with a map
	std::map<bp_domain, long double> probabilities;
	/// Set default value at initialization time
	long double nullval;
public:
	// Constructor
	/// Construct a new message with all elements initialized to nullval
	bp_message(long double nullval);
	// Set & access mapped values
	long double& operator[](const bp_domain& value);
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
	void normalize(int domain_size);
	// Extract maximum value
	bp_domain extract_max();
};

#endif
