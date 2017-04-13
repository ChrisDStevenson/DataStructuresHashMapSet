// Submitter: cdsteve1(Stevenson, Christopher)
// Partner: kpmckeow(Mckeown, Kelly)
// We certify that we worked cooperatively on this programming
//   assignment, according to the rules for pair programming
#ifndef HASH_SET_HPP_
#define HASH_SET_HPP_

#include <string>
#include <iostream>
#include <sstream>
#include <initializer_list>
#include "ics_exceptions.hpp"
#include "pair.hpp"


namespace ics {


#ifndef undefinedhashdefined
#define undefinedhashdefined
template<class T>
int undefinedhash (const T& a) {return 0;}
#endif /* undefinedhashdefined */

//Instantiate the templated class supplying thash(a): produces a hash value for a.
//If thash is defaulted to undefinedhash in the template, then a constructor must supply chash.
//If both thash and chash are supplied, then they must be the same (by ==) function.
//If neither is supplied, or both are supplied but different, TemplateFunctionError is raised.
//The (unique) non-undefinedhash value supplied by thash/chash is stored in the instance variable hash.
template<class T, int (*thash)(const T& a) = undefinedhash<T>> class HashSet {
  public:
    typedef int (*hashfunc) (const T& a);

    //Destructor/Constructors
    ~HashSet ();

    HashSet (double the_load_threshold = 1.0, int (*chash)(const T& a) = undefinedhash<T>);
    explicit HashSet (int initial_bins, double the_load_threshold = 1.0, int (*chash)(const T& k) = undefinedhash<T>);
    HashSet (const HashSet<T,thash>& to_copy, double the_load_threshold = 1.0, int (*chash)(const T& a) = undefinedhash<T>);
    explicit HashSet (const std::initializer_list<T>& il, double the_load_threshold = 1.0, int (*chash)(const T& a) = undefinedhash<T>);

    //Iterable class must support "for-each" loop: .begin()/.end() and prefix ++ on returned result
    template <class Iterable>
    explicit HashSet (const Iterable& i, double the_load_threshold = 1.0, int (*chash)(const T& a) = undefinedhash<T>);


    //Queries
    bool empty      () const;
    int  size       () const;
    bool contains   (const T& element) const;
    std::string str () const; //supplies useful debugging information; contrast to operator <<

    //Iterable class must support "for-each" loop: .begin()/.end() and prefix ++ on returned result
    template <class Iterable>
    bool contains_all (const Iterable& i) const;


    //Commands
    int  insert (const T& element);
    int  erase  (const T& element);
    void clear  ();

    //Iterable class must support "for" loop: .begin()/.end() and prefix ++ on returned result

    template <class Iterable>
    int insert_all(const Iterable& i);

    template <class Iterable>
    int erase_all(const Iterable& i);

    template<class Iterable>
    int retain_all(const Iterable& i);


    //Operators
    HashSet<T,thash>& operator = (const HashSet<T,thash>& rhs);
    bool operator == (const HashSet<T,thash>& rhs) const;
    bool operator != (const HashSet<T,thash>& rhs) const;
    bool operator <= (const HashSet<T,thash>& rhs) const;
    bool operator <  (const HashSet<T,thash>& rhs) const;
    bool operator >= (const HashSet<T,thash>& rhs) const;
    bool operator >  (const HashSet<T,thash>& rhs) const;

    template<class T2, int (*hash2)(const T2& a)>
    friend std::ostream& operator << (std::ostream& outs, const HashSet<T2,hash2>& s);



  private:
    class LN;

  public:
    class Iterator {
      public:
        typedef pair<int,LN*> Cursor;

        //Private constructor called in begin/end, which are friends of HashSet<T,thash>
        ~Iterator();
        T           erase();
        std::string str  () const;
        HashSet<T,thash>::Iterator& operator ++ ();
        HashSet<T,thash>::Iterator  operator ++ (int);
        bool operator == (const HashSet<T,thash>::Iterator& rhs) const;
        bool operator != (const HashSet<T,thash>::Iterator& rhs) const;
        T& operator *  () const;
        T* operator -> () const;
        friend std::ostream& operator << (std::ostream& outs, const HashSet<T,thash>::Iterator& i) {
          outs << i.str(); //Use the same meaning as the debugging .str() method
          return outs;
        }
        friend Iterator HashSet<T,thash>::begin () const;
        friend Iterator HashSet<T,thash>::end   () const;

      private:
        //If can_erase is false, current indexes the "next" value (must ++ to reach it)
        Cursor              current; //Bin Index + LN* pointer; stops if LN* == nullptr
        HashSet<T,thash>*   ref_set;
        int                 expected_mod_count;
        bool                can_erase = true;

        //Helper methods
        void advance_cursors();

        //Called in friends begin/end
        Iterator(HashSet<T,thash>* iterate_over, bool from_begin);
    };


    Iterator begin () const;
    Iterator end   () const;


  private:
    class LN {
      public:
        LN ()                      {}
        LN (const LN& ln)          : value(ln.value), next(ln.next){}
        LN (T v,  LN* n = nullptr) : value(v), next(n){}

        T   value;
        LN* next   = nullptr;
    };

public:
  int (*hash)(const T& k);   //Hashing function used (from template or constructor)
private:
  LN** set      = nullptr;   //Pointer to array of pointers: each bin stores a list with a trailer node
  double load_threshold;     //used/bins <= load_threshold
  int bins      = 1;         //# bins in array (should start >= 1 so hash_compress doesn't divide by 0)
  int used      = 0;         //Cache for number of key->value pairs in the hash table
  int mod_count = 0;         //For sensing concurrent modification


  //Helper methods
  int   hash_compress        (const T& key)              const;  //hash function ranged to [0,bins-1]
  LN*   find_element         (const T& element)          const;  //Returns reference to element's node or nullptr
  bool  has_value            (LN* l, T value) const;

  LN*   copy_list            (LN*   l)                   const;  //Copy the elements in a bin (order irrelevant)
  LN**  copy_hash_table      (LN** ht, int bins)         const;  //Copy the bins/keys/values in ht tree (order in bins irrelevant)
  void  rehash_all          (LN*& iter);                       //Rehash/move all  keys/values


  void  ensure_load_threshold(int new_used);                     //Reallocate if load_threshold > load_threshold
  void  delete_hash_table    (LN**& ht, int bins);               //Deallocate all LN in ht (and the ht itself; ht == nullptr)
  LN**  new_hash_table (int bins);                             //Create a new hash table with n bins

    };





//HashSet class and related definitions

////////////////////////////////////////////////////////////////////////////////
//
//Destructor/Constructors

template<class T, int (*thash)(const T& a)>
HashSet<T,thash>::~HashSet() {
}


template<class T, int (*thash)(const T& a)>
HashSet<T,thash>::HashSet(double the_load_threshold, int (*chash)(const T& element))
        : hash(thash != (hashfunc)undefinedhash<T> ? thash: chash), load_threshold(the_load_threshold){
    if (hash == (hashfunc)undefinedhash<T>)
        throw TemplateFunctionError("HashSet::default constructor: neither specified");
    if (thash != (hashfunc)undefinedhash<T> && chash != (hashfunc)undefinedhash<T> && thash != chash)
        throw TemplateFunctionError("HashSet::default constructor: both specified and different");
    set = new_hash_table(1);
}


template<class T, int (*thash)(const T& a)>
HashSet<T,thash>::HashSet(int initial_bins, double the_load_threshold, int (*chash)(const T& element))
        : hash(thash != (hashfunc)undefinedhash<T> ? thash: chash), load_threshold(the_load_threshold), bins(initial_bins) {
    if (hash == (hashfunc)undefinedhash<T>)
        throw TemplateFunctionError("HashSet::default constructor: neither specified");
    if (thash != (hashfunc)undefinedhash<T> && chash != (hashfunc)undefinedhash<T> && thash != chash)
        throw TemplateFunctionError("HashSet::default constructor: both specified and different");

    set = new_hash_table(initial_bins);
}


template<class T, int (*thash)(const T& a)>
HashSet<T,thash>::HashSet(const HashSet<T,thash>& to_copy, double the_load_threshold, int (*chash)(const T& element))
        : hash(thash != (hashfunc)undefinedhash<T> ? thash: chash), load_threshold(the_load_threshold){

    if (hash == (hashfunc)undefinedhash<T>)
        hash = to_copy.hash;
    if (thash != (hashfunc)undefinedhash<T> && chash != (hashfunc)undefinedhash<T> && thash != chash)
        throw TemplateFunctionError("HashSet::copy constructor: both specified and different");

    if (thash != chash){
        hash = to_copy.hash;
        set = new_hash_table(1);
        auto iter = to_copy.begin();
        while ( iter != to_copy.end()) {
            insert(*iter);
            iter++;
        }
    }
    else{
        set = copy_hash_table(to_copy.set, to_copy.bins);
        bins = to_copy.bins;
    }
    used = to_copy.used;
}


template<class T, int (*thash)(const T& a)>
HashSet<T,thash>::HashSet(const std::initializer_list<T>& il, double the_load_threshold, int (*chash)(const T& element))
        : hash(thash != (hashfunc)undefinedhash<T> ? thash: chash), load_threshold(the_load_threshold){
    set = new_hash_table(1);
    for (const T& pq_elem : il)
        insert(pq_elem);
}


template<class T, int (*thash)(const T& a)>
template<class Iterable>
HashSet<T,thash>::HashSet(const Iterable& i, double the_load_threshold, int (*chash)(const T& a))
        : hash(thash != (hashfunc)undefinedhash<T> ? thash: chash), load_threshold(the_load_threshold){
    set = new_hash_table(1);
    for (const T& pq_elem : i)
        insert(pq_elem);
}


////////////////////////////////////////////////////////////////////////////////
//
//Queries

template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::empty() const {
    return used == 0;
}


template<class T, int (*thash)(const T& a)>
int HashSet<T,thash>::size() const {
    return used;
}


template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::contains (const T& element) const {
    if (!empty())
        return has_value(set[hash_compress(element)], element);
    return false;
}


template<class T, int (*thash)(const T& a)>
std::string HashSet<T,thash>::str() const {
    std::ostringstream answer;
    for (int i = 0; i < bins; ++i) {
        answer << "bin[" << i << "]:";
        for (LN *iter = set[i]; iter != nullptr; iter = iter->next) {
            if (iter->next != nullptr)
                answer << " " << iter->value;
            else
                answer << " TRAILER";
        }
        answer << std::endl;
    }
    return answer.str();
}


template<class T, int (*thash)(const T& a)>
template <class Iterable>
bool HashSet<T,thash>::contains_all(const Iterable& i) const {
    for (const T& v : i){
        if (!contains(v))
            return false;
    }
    return true;
}


////////////////////////////////////////////////////////////////////////////////
//
//Commands

template<class T, int (*thash)(const T& a)>
int HashSet<T,thash>::insert(const T& element) {
    LN * key_exists = (!empty()) ? find_element(element) : nullptr;
    if (key_exists != nullptr){
        return 0;
    }
    else {
        ensure_load_threshold(used + 1);
        int location = hash_compress(element);
        set[location] = new LN(element, set[location]);
        ++used;
    }
    ++mod_count;
    return 1;
}


template<class T, int (*thash)(const T& a)>
int HashSet<T,thash>::erase(const T& element) {

    if (empty())
        return 0;

    LN* to_erase = find_element(element);

    if (to_erase == nullptr)
        return 0;

    T value = to_erase->value;
    int index = hash_compress(element);

    LN* current = set[index];

    if (current == to_erase){
        set[index] = current->next;
    }
    else{
        while (current->next != to_erase)
            current = current->next;
        current->next = to_erase->next;
    }
    delete to_erase;
    --used;
    ++mod_count;
    return 1;
}


template<class T, int (*thash)(const T& a)>
void HashSet<T,thash>::clear() {
    used = 0;
    ++mod_count;
    delete_hash_table(set,bins);
    bins = 1;
    load_threshold = 0;
    set = new_hash_table(1);
}


template<class T, int (*thash)(const T& a)>
template<class Iterable>
int HashSet<T,thash>::insert_all(const Iterable& i) {
    int count = 0;
    for (const T & iter : i){
        insert(iter);
        count++;
    }
    return count;
}


template<class T, int (*thash)(const T& a)>
template<class Iterable>
int HashSet<T,thash>::erase_all(const Iterable& i) {
    int count = 0;
    for (const T & iter : i){
        erase(iter);
        count++;
    }
    return count;
}


template<class T, int (*thash)(const T& a)>
template<class Iterable>
int HashSet<T,thash>::retain_all(const Iterable& i) {
    HashSet s(i);
    int count = 0;
    
    for (auto iter = begin(); iter != end(); iter++)
        if (!s.contains(*iter)) {
            iter.erase();
            ++count;
        }

    return count;
}


////////////////////////////////////////////////////////////////////////////////
//
//Operators

template<class T, int (*thash)(const T& a)>
HashSet<T,thash>& HashSet<T,thash>::operator = (const HashSet<T,thash>& rhs) {
    clear();
    hash = rhs.hash;
    set = copy_hash_table(rhs.set, rhs.bins);
    used = rhs.size();
    ++mod_count;
    return *this;
}


template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::operator == (const HashSet<T,thash>& rhs) const {
    if (used == rhs.size()) {
        if (this == &rhs)
            return true;
        if (empty() && rhs.empty())
            return true;
        for (auto iter = begin(); iter != end(); ++iter) {
            LN *find_relement = rhs.find_element(*iter);
            if (find_relement == nullptr || find_relement->value != *iter)
                return false;
        }
        return true;
    }
    return false;
}


template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::operator != (const HashSet<T,thash>& rhs) const {
    return !(*this == rhs);
}


template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::operator <= (const HashSet<T,thash>& rhs) const {
    if (this == &rhs)
        return true;

    if (used > rhs.size())
        return false;

    for (auto iter = begin(); iter != end(); ++iter)
        if (!rhs.contains(*iter))
            return false;

    return true;
}

template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::operator < (const HashSet<T,thash>& rhs) const {
    if (this == &rhs)
        return false;

    if (used >= rhs.size())
        return false;

    for (auto iter = begin(); iter != end(); ++iter)
        if (!rhs.contains(*iter))
            return false;

    return true;
}


template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::operator >= (const HashSet<T,thash>& rhs) const {
    return rhs <= *this;
}


template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::operator > (const HashSet<T,thash>& rhs) const {
    return rhs < *this;
}


template<class T, int (*thash)(const T& a)>
std::ostream& operator << (std::ostream& outs, const HashSet<T,thash>& s) {
    outs << "set[";
    auto iter = s.begin();
    while (iter != s.end()){
        outs << *iter;
        iter++;
        if (iter != s.end())
            outs << ", ";
    }
    outs << "]";
    return outs;
}


////////////////////////////////////////////////////////////////////////////////
//
//Iterator constructors

template<class T, int (*thash)(const T& a)>
auto HashSet<T,thash>::begin () const -> HashSet<T,thash>::Iterator {
    return Iterator(const_cast<HashSet<T,thash>*>(this),true);
}


template<class T, int (*thash)(const T& a)>
auto HashSet<T,thash>::end () const -> HashSet<T,thash>::Iterator {
    return Iterator(const_cast<HashSet<T,thash>*>(this),false);
}


///////////////////////////////////////////////////////////////////////////////
//
//Private helper methods

template<class T, int (*thash)(const T& a)>
int HashSet<T,thash>::hash_compress (const T& element) const {
    return std::abs(hash(element)) % bins;
}


template<class T, int (*thash)(const T& a)>
typename HashSet<T,thash>::LN* HashSet<T,thash>::find_element (const T& element) const {
    int index = hash_compress(element);
    LN* iter = set[hash_compress(element)];
    while (iter->next != nullptr){
        if (iter->value == element)
            return iter;
        iter = iter->next;
    }
    return nullptr;
}


template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::has_value (LN* l, T value) const{
    if (l->next == nullptr)
        return false;
    if (l->value == value)
        return true;
    return has_value(l->next, value);
};


    template<class T, int (*thash)(const T& a)>
typename HashSet<T,thash>::LN* HashSet<T,thash>::copy_list (LN* l) const {
    if (l == nullptr)
        return nullptr;
    return new LN(l->value, copy_list(l->next));
}


template<class T, int (*thash)(const T& a)>
typename HashSet<T,thash>::LN** HashSet<T,thash>::copy_hash_table (LN** ht, int bins) const {
    LN** copy_h = new LN*[bins];
    for (int i = 0; i < bins; ++i)
        copy_h[i] = copy_list(ht[i]);
    return copy_h;
}


template<class T, int (*thash)(const T& a)>
void HashSet<T,thash>::ensure_load_threshold(int new_used) {
    double load_threshold = new_used / bins;
    if (load_threshold > 1){
        LN ** old_h = set;
        set = new_hash_table(bins * 2);

        int old_bins = bins;
        bins = bins * 2;

        for (int i = 0; i < old_bins; ++i)
        {
            rehash_all(old_h[i]);
        }
        delete[] old_h;
    }
}

template<class T, int (*thash)(const T& a)>
void HashSet<T,thash>::rehash_all (LN* & iter) {
    if (iter->next == nullptr) {
        delete iter; //Remove trailer after reassigning all values
        return;
    }
    int dest_bin = hash_compress(iter->value);
    LN* next = iter->next;
    iter->next = set[dest_bin];
    set[dest_bin] = iter;
    rehash_all(next);
}



template<class T, int (*thash)(const T& a)>
typename HashSet<T,thash>:: LN** HashSet<T,thash>::new_hash_table (int bins) {
    LN ** new_h = new LN*[bins];
    for (int i = 0; i < bins; ++i)
        new_h[i] = new LN();
    return new_h;
}
    
template<class T, int (*thash)(const T& a)>
void HashSet<T,thash>::delete_hash_table (LN**& ht, int bins) {
    if (!empty()) {
        for (int i = 0; i < bins; ++i) {
            LN *iter = ht[i];
            while (iter != nullptr) {
                LN *toDelete = iter;
                iter = iter->next;
                delete toDelete;
            }
        }
    }
    delete[] ht;
}






////////////////////////////////////////////////////////////////////////////////
//
//Iterator class definitions

template<class T, int (*thash)(const T& a)>
void HashSet<T,thash>::Iterator::advance_cursors() {
    while (true){
        if (current.first == ref_set->bins - 1 && current.second->next == nullptr){
            current = Cursor(-1, nullptr);
            return;
        }
        if (current.second->next == nullptr){
            current.second = ref_set->set[++current.first];
            if (current.second->next != nullptr)
                return;
        }
        else{
            current.second = current.second->next;
            if (current.second->next != nullptr)
                return;
        }
    }
}


template<class T, int (*thash)(const T& a)>
HashSet<T,thash>::Iterator::Iterator(HashSet<T,thash>* iterate_over, bool begin)
        : ref_set(iterate_over), expected_mod_count(ref_set->mod_count) {

    current = (begin) ? Cursor(0, ref_set->set[0]) : Cursor(-1, nullptr);
    if (begin && current.second->next == nullptr)
        advance_cursors();

}


template<class T, int (*thash)(const T& a)>
HashSet<T,thash>::Iterator::~Iterator()
{}


template<class T, int (*thash)(const T& a)>
T HashSet<T,thash>::Iterator::erase() {
    if (expected_mod_count != ref_set->mod_count)
        throw ConcurrentModificationError("HashSet::Iterator::erase");
    if (!can_erase)
        throw CannotEraseError("HashSet::Iterator::erase Iterator cursor already erased");
    if (current.first == -1 && current.second == nullptr)
        throw CannotEraseError("HashSet::Iterator::erase Iterator cursor beyond data structure");

    can_erase = false;
    T to_return = current.second->value;
    advance_cursors();
    ref_set->erase(to_return);
    expected_mod_count = ref_set->mod_count;
    return to_return;
}


template<class T, int (*thash)(const T& a)>
std::string HashSet<T,thash>::Iterator::str() const {
    std::ostringstream answer;
    answer << ref_set->str() << "/current=" << current.first  << "/expected_mod_count=" << expected_mod_count << "/can_erase=" << can_erase;
    return answer.str();
}


template<class T, int (*thash)(const T& a)>
auto  HashSet<T,thash>::Iterator::operator ++ () -> HashSet<T,thash>::Iterator& {
    if (expected_mod_count != ref_set->mod_count)
        throw ConcurrentModificationError("HashSet::Iterator::operator ++(int)");

    if (current.first == -1 && current.second == nullptr)
        return *this;

    if (can_erase)
        advance_cursors();
    else
        can_erase = true;

    return *this;
}


template<class T, int (*thash)(const T& a)>
auto  HashSet<T,thash>::Iterator::operator ++ (int) -> HashSet<T,thash>::Iterator {
    if (expected_mod_count != ref_set->mod_count)
        throw ConcurrentModificationError("HashSet::Iterator::operator ++(int)");

    if (current.first == -1 && current.second == nullptr)
        return *this;

    Iterator to_return(*this);
    if (can_erase)
        advance_cursors();
    else
        can_erase = true;

    return to_return;
}


template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::Iterator::operator == (const HashSet<T,thash>::Iterator& rhs) const {
    const Iterator* rhsASI = dynamic_cast<const Iterator*>(&rhs);
    if (rhsASI == 0)
        throw IteratorTypeError("HashSet::Iterator::operator ==");
    if (expected_mod_count != ref_set->mod_count)
        throw ConcurrentModificationError("HashSet::Iterator::operator ==");
    if (ref_set != rhsASI->ref_set)
        throw ComparingDifferentIteratorsError("HashSet::Iterator::operator ==");
    return (current == rhsASI->current);
}


template<class T, int (*thash)(const T& a)>
bool HashSet<T,thash>::Iterator::operator != (const HashSet<T,thash>::Iterator& rhs) const {
    return !(*this == rhs);
}

template<class T, int (*thash)(const T& a)>
T& HashSet<T,thash>::Iterator::operator *() const {
    if (expected_mod_count != ref_set->mod_count)
        throw ConcurrentModificationError("HashSet::Iterator::operator *");
    if (!can_erase || current == Cursor(-1, nullptr)) {
        std::ostringstream where;
        where << " when size = " << ref_set->size();
        throw IteratorPositionIllegal("HashSet::Iterator::operator * Iterator illegal: "+where.str());
    }
    return current.second->value;
}

template<class T, int (*thash)(const T& a)>
T* HashSet<T,thash>::Iterator::operator ->() const {
    if (expected_mod_count != ref_set->mod_count)
        throw ConcurrentModificationError("HashSet::Iterator::operator ->");
    if (!can_erase || ref_set->used == 0){
        std::ostringstream where;
        where << current.second->value << " when size = " << ref_set->size();
        throw IteratorPositionIllegal("HashSet::Iterator::operator -> Iterator illegal: "+where.str());
    }
    return &current.second->value;
}

}

#endif /* HASH_SET_HPP_ */
