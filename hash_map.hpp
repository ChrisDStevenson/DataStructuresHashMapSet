// Submitter: cdsteve1(Stevenson, Christopher)
// Partner: kpmckeow(Mckeown, Kelly)
// We certify that we worked cooperatively on this programming
//   assignment, according to the rules for pair programming
#ifndef HASH_MAP_HPP_
#define HASH_MAP_HPP_

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
template<class KEY,class T, int (*thash)(const KEY& a) = undefinedhash<KEY>> class HashMap {
  public:
    typedef ics::pair<KEY,T>   Entry;
    typedef int (*hashfunc) (const KEY& a);

    //Destructor/Constructors
    ~HashMap ();

    HashMap          (double the_load_threshold = 1.0, int (*chash)(const KEY& a) = undefinedhash<KEY>);
    explicit HashMap (int initial_bins, double the_load_threshold = 1.0, int (*chash)(const KEY& k) = undefinedhash<KEY>);
    HashMap          (const HashMap<KEY,T,thash>& to_copy, double the_load_threshold = 1.0, int (*chash)(const KEY& a) = undefinedhash<KEY>);
    explicit HashMap (const std::initializer_list<Entry>& il, double the_load_threshold = 1.0, int (*chash)(const KEY& a) = undefinedhash<KEY>);

    //Iterable class must support "for-each" loop: .begin()/.end() and prefix ++ on returned result
    template <class Iterable>
    explicit HashMap (const Iterable& i, double the_load_threshold = 1.0, int (*chash)(const KEY& a) = undefinedhash<KEY>);


    //Queries
    bool empty      () const;
    int  size       () const;
    bool has_key    (const KEY& key) const;
    bool has_value  (const T& value) const;
    std::string str () const; //supplies useful debugging information; contrast to operator <<


    //Commands
    T    put   (const KEY& key, const T& value);
    T    erase (const KEY& key);
    void clear ();

    //Iterable class must support "for-each" loop: .begin()/.end() and prefix ++ on returned result
    template <class Iterable>
    int put_all(const Iterable& i);


    //Operators

    T&       operator [] (const KEY&);
    const T& operator [] (const KEY&) const;
    HashMap<KEY,T,thash>& operator = (const HashMap<KEY,T,thash>& rhs);
    bool operator == (const HashMap<KEY,T,thash>& rhs) const;
    bool operator != (const HashMap<KEY,T,thash>& rhs) const;

    template<class KEY2,class T2, int (*hash2)(const KEY2& a)>
    friend std::ostream& operator << (std::ostream& outs, const HashMap<KEY2,T2,hash2>& m);



  private:
    class LN;

  public:
    class Iterator {
      public:
         typedef pair<int,LN*> Cursor;

        //Private constructor called in begin/end, which are friends of HashMap<T>
        ~Iterator();
        Entry       erase();
        std::string str  () const;
        HashMap<KEY,T,thash>::Iterator& operator ++ ();
        HashMap<KEY,T,thash>::Iterator  operator ++ (int);
        bool operator == (const HashMap<KEY,T,thash>::Iterator& rhs) const;
        bool operator != (const HashMap<KEY,T,thash>::Iterator& rhs) const;
        Entry& operator *  () const;
        Entry* operator -> () const;
        friend std::ostream& operator << (std::ostream& outs, const HashMap<KEY,T,thash>::Iterator& i) {
          outs << i.str(); //Use the same meaning as the debugging .str() method
          return outs;
        }
        friend Iterator HashMap<KEY,T,thash>::begin () const;
        friend Iterator HashMap<KEY,T,thash>::end   () const;

      private:
        //If can_erase is false, current indexes the "next" value (must ++ to reach it)
        Cursor                current; //Bin Index + LN* pointer; stops if LN* == nullptr
        HashMap<KEY,T,thash>* ref_map;
        int                   expected_mod_count;
        bool                  can_erase = true;

        //Helper methods
        void advance_cursors();

        //Called in friends begin/end
        Iterator(HashMap<KEY,T,thash>* iterate_over, bool from_begin);
    };


    Iterator begin () const;
    Iterator end   () const;


  private:
    class LN {
    public:
      LN ()                         : next(nullptr){}
      LN (const LN& ln)             : value(ln.value), next(ln.next){}
      LN (Entry v, LN* n = nullptr) : value(v), next(n){}

      Entry value;
      LN*   next;
  };

  int (*hash)(const KEY& k);  //Hashing function used (from template or constructor)
  LN** map      = nullptr;    //Pointer to array of pointers: each bin stores a list with a trailer node
  double load_threshold;      //used/bins <= load_threshold
  int bins      = 1;          //# bins in array (should start >= 1 so hash_compress doesn't divide by 0)
  int used      = 0;          //Cache for number of key->value pairs in the hash table
  int mod_count = 0;          //For sensing concurrent modification


  //Helper methods
  int   hash_compress        (const KEY& key)          const;  //hash function ranged to [0,bins-1]
  LN*   find_key             (const KEY& key)          const;  //Returns reference to key's node or nullptr
  bool  has_value            (LN* l, T value) const;

  LN*   copy_list            (LN*   l)                 const;  //Copy the keys/values in a bin (order irrelevant)
  LN**  copy_hash_table      (LN** ht, int bins)       const;  //Copy the bins/keys/values in ht tree (order in bins irrelevant
  void  rehash_all          (LN*& iter);                       //Rehash/move all  keys/values


  void  ensure_load_threshold(int new_used);                   //Reallocate if load_factor > load_threshold
  void  delete_hash_table    (LN**& ht, int bins);             //Deallocate all LN in ht (and the ht itself; ht == nullptr)
  LN**  new_hash_table (int bins);                             //Create a new hash table with n bins

};





////////////////////////////////////////////////////////////////////////////////
//
//HashMap class and related definitions

//Destructor/Constructors

template<class KEY,class T, int (*thash)(const KEY& a)>
HashMap<KEY,T,thash>::~HashMap() {
    delete_hash_table(map, bins);
}


template<class KEY,class T, int (*thash)(const KEY& a)>
HashMap<KEY,T,thash>::HashMap(double the_load_threshold, int (*chash)(const KEY& k))
: hash(thash != (hashfunc)undefinedhash<KEY> ? thash: chash), load_threshold(the_load_threshold){
    if (hash == (hashfunc)undefinedhash<KEY>)
        throw TemplateFunctionError("HashMap::default constructor: neither specified");
    if (thash != (hashfunc)undefinedhash<KEY> && chash != (hashfunc)undefinedhash<KEY> && thash != chash)
        throw TemplateFunctionError("HashMap::default constructor: both specified and different");
    map = new_hash_table(1);

}


template<class KEY,class T, int (*thash)(const KEY& a)>
HashMap<KEY,T,thash>::HashMap(int initial_bins, double the_load_threshold, int (*chash)(const KEY& k))
        : hash(thash != (hashfunc)undefinedhash<KEY> ? thash: chash), load_threshold(the_load_threshold), bins(initial_bins) {
    if (hash == (hashfunc)undefinedhash<KEY>)
        throw TemplateFunctionError("HashMap::default constructor: neither specified");
    if (thash != (hashfunc)undefinedhash<KEY> && chash != (hashfunc)undefinedhash<KEY> && thash != chash)
        throw TemplateFunctionError("HashMap::default constructor: both specified and different");

    map = new_hash_table(initial_bins);
}


template<class KEY,class T, int (*thash)(const KEY& a)>
HashMap<KEY,T,thash>::HashMap(const HashMap<KEY,T,thash>& to_copy, double the_load_threshold, int (*chash)(const KEY& a))
 : hash(thash != (hashfunc)undefinedhash<KEY> ? thash: chash), load_threshold(the_load_threshold){

    if (hash == (hashfunc)undefinedhash<KEY>)
        hash = to_copy.hash;
    if (thash != (hashfunc)undefinedhash<KEY> && chash != (hashfunc)undefinedhash<KEY> && thash != chash)
        throw TemplateFunctionError("HashMap::copy constructor: both specified and different");

    if (thash != chash){
        hash = to_copy.hash;
        map = new_hash_table(1);
        auto iter = to_copy.begin();
        while ( iter != to_copy.end()) {
            put(iter->first, iter->second);
            iter++;
        }
    }
    else{
        map = copy_hash_table(to_copy.map, to_copy.bins);
        bins = to_copy.bins;
    }
    used = to_copy.used;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
HashMap<KEY,T,thash>::HashMap(const std::initializer_list<Entry>& il, double the_load_threshold, int (*chash)(const KEY& k))
        : hash(thash != (hashfunc)undefinedhash<KEY> ? thash: chash), load_threshold(the_load_threshold){
    map = new_hash_table(1);
    for (const Entry& pq_elem : il)
        put(pq_elem.first, pq_elem.second);
}


template<class KEY,class T, int (*thash)(const KEY& a)>
template <class Iterable>
HashMap<KEY,T,thash>::HashMap(const Iterable& i, double the_load_threshold, int (*chash)(const KEY& k))
        : hash(thash != (hashfunc)undefinedhash<KEY> ? thash: chash), load_threshold(the_load_threshold){
    map = new_hash_table(1);
    for (const Entry& pq_elem : i)
        put(pq_elem.first, pq_elem.second);
}


////////////////////////////////////////////////////////////////////////////////
//
//Queries

template<class KEY,class T, int (*thash)(const KEY& a)>
bool HashMap<KEY,T,thash>::empty() const {
    return used == 0;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
int HashMap<KEY,T,thash>::size() const {
    return used;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
bool HashMap<KEY,T,thash>::has_key (const KEY& key) const {
    return (!empty() && find_key(key) != nullptr);
}


template<class KEY,class T, int (*thash)(const KEY& a)>
bool HashMap<KEY,T,thash>::has_value (const T& value) const {
    if (!empty())
        for (int i = 0; i < bins; ++i)
            if (has_value(map[i], value))
                return true;
    return false;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
std::string HashMap<KEY,T,thash>::str() const {
    std::ostringstream answer;
    for (int i = 0; i < bins; ++i) {
        answer << "bin[" << i << "]:";
        for (LN *iter = map[i]; iter != nullptr; iter = iter->next) {
            if (iter->next != nullptr)
                answer << " " << iter->value.first << "->" << iter->value.second;
            else
                answer << " TRAILER";
        }
        answer << std::endl;
    }
    return answer.str();

}


////////////////////////////////////////////////////////////////////////////////
//
//Commands

template<class KEY,class T, int (*thash)(const KEY& a)>
T HashMap<KEY,T,thash>::put(const KEY& key, const T& value) {
    T return_value = value;
    LN * key_exists = (!empty()) ? find_key(key) : nullptr;
    if (key_exists != nullptr){
        return_value = key_exists->value.second;
        key_exists->value.second = value;
    }
    else {
        ensure_load_threshold(used + 1);
        int location = hash_compress(key);
        map[location] = new LN(Entry(key, value), map[location]);
        ++used;
    }
    ++mod_count;
    return return_value;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
T HashMap<KEY,T,thash>::erase(const KEY& key) {
    if (empty())
        throw ics::KeyError("Specified key was not found");

    LN* to_erase = find_key(key);

    if (to_erase == nullptr)
        throw ics::KeyError("Specified key was not found");

    T value = to_erase->value.second;
    int index = hash_compress(key);

    LN* current = map[index];

    if (current == to_erase){
        map[index] = current->next;
    }
    else{
        while (current->next != to_erase)
            current = current->next;
        current->next = to_erase->next;
    }
    delete to_erase;
    --used;
    ++mod_count;
    return value;

}


template<class KEY,class T, int (*thash)(const KEY& a)>
void HashMap<KEY,T,thash>::clear() {
    used = 0;
    ++mod_count;
    delete_hash_table(map,bins);
    bins = 1;
    load_threshold = 0;
    map = new_hash_table(1);
}


template<class KEY,class T, int (*thash)(const KEY& a)>
template<class Iterable>
int HashMap<KEY,T,thash>::put_all(const Iterable& i) {
    int count = 0;
    for (const Entry & iter : i){
        put(iter.first,iter.second);
        count++;
    }
    return count;
}


////////////////////////////////////////////////////////////////////////////////
//
//Operators

template<class KEY,class T, int (*thash)(const KEY& a)>
T& HashMap<KEY,T,thash>::operator [] (const KEY& key) {
    LN * search_map = nullptr;
    if (!empty())
        search_map = find_key(key);
    if (search_map == nullptr){
        ensure_load_threshold(used + 1);
        int location = hash_compress(key);
        LN * NewItem = new LN(Entry(key, T()), map[location]);
        map[location] = NewItem;
        ++used;
        ++mod_count;
        return NewItem->value.second;
    }
    else
        return search_map->value.second;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
const T& HashMap<KEY,T,thash>::operator [] (const KEY& key) const {
    const LN * find_pair = find_key(key);
    return find_pair->value.second;

}


template<class KEY,class T, int (*thash)(const KEY& a)>
HashMap<KEY,T,thash>& HashMap<KEY,T,thash>::operator = (const HashMap<KEY,T,thash>& rhs) {
    clear();
    hash = rhs.hash;
    map = copy_hash_table(rhs.map, rhs.bins);
    used = rhs.size();
    ++mod_count;
    return *this;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
bool HashMap<KEY,T,thash>::operator == (const HashMap<KEY,T,thash>& rhs) const {
    if (used == rhs.size()) {
        if (this == &rhs)
            return true;
        if (empty() && rhs.empty())
            return true;
        for (auto iter = begin(); iter != end(); ++iter) {
            LN *find_rkey = rhs.find_key(iter->first);
            if (find_rkey == nullptr || find_rkey->value.second != iter->second)
                return false;
        }
        return true;
    }
    return false;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
bool HashMap<KEY,T,thash>::operator != (const HashMap<KEY,T,thash>& rhs) const {
    return !(*this == rhs);
}


template<class KEY,class T, int (*thash)(const KEY& a)>
std::ostream& operator << (std::ostream& outs, const HashMap<KEY,T,thash>& m) {
    outs << "map[";
    auto iter = m.begin();
    while (iter != m.end()){
        outs << iter->first << "->" << iter->second;
        iter++;
        if (iter != m.end())
            outs << ", ";
    }
    outs << "]";
    return outs;
}


////////////////////////////////////////////////////////////////////////////////
//
//Iterator constructors

template<class KEY,class T, int (*thash)(const KEY& a)>
auto HashMap<KEY,T,thash>::begin () const -> HashMap<KEY,T,thash>::Iterator {
    return Iterator(const_cast<HashMap<KEY,T,thash>*>(this),true);
}


template<class KEY,class T, int (*thash)(const KEY& a)>
auto HashMap<KEY,T,thash>::end () const -> HashMap<KEY,T,thash>::Iterator {
    return Iterator(const_cast<HashMap<KEY,T,thash>*>(this),false);
}


///////////////////////////////////////////////////////////////////////////////
//
//Private helper methods

template<class KEY,class T, int (*thash)(const KEY& a)>
int HashMap<KEY,T,thash>::hash_compress (const KEY& key) const {
    return std::abs(hash(key)) % bins;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
typename HashMap<KEY,T,thash>::LN* HashMap<KEY,T,thash>::find_key (const KEY& key) const{
    int index = hash_compress(key);
    LN* iter = map[hash_compress(key)];
    while (iter->next != nullptr){
        if (iter->value.first == key)
            return iter;
        iter = iter->next;
    }
    return nullptr;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
bool HashMap<KEY,T,thash>::has_value (LN* l, T value) const{
    if (l->next == nullptr)
        return false;
    if (l->value.second == value)
        return true;
    return has_value(l->next, value);
};

template<class KEY,class T, int (*thash)(const KEY& a)>
typename HashMap<KEY,T,thash>::LN* HashMap<KEY,T,thash>::copy_list (LN* l) const {
    if (l == nullptr)
        return nullptr;
    return new LN(l->value, copy_list(l->next));
}


template<class KEY,class T, int (*thash)(const KEY& a)>
typename HashMap<KEY,T,thash>::LN** HashMap<KEY,T,thash>::copy_hash_table (LN** ht, int bins) const {
    LN** copy_h = new LN*[bins];
    for (int i = 0; i < bins; ++i)
        copy_h[i] = copy_list(ht[i]);
    return copy_h;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
void HashMap<KEY,T,thash>::ensure_load_threshold(int new_used) {
    double load_threshold = new_used / bins;
    if (load_threshold > 1){
        LN ** old_h = map;
        map = new_hash_table(bins * 2);

        int old_bins = bins;
        bins = bins * 2;

        for (int i = 0; i < old_bins; ++i)
        {
            rehash_all(old_h[i]);
        }
        delete[] old_h;
    }
}


template<class KEY,class T, int (*thash)(const KEY& a)>
void HashMap<KEY,T,thash>::rehash_all (LN* & iter) {
    if (iter->next == nullptr) {
        delete iter; //Remove trailer after reassigning all values
        return;
    }
    int dest_bin = hash_compress(iter->value.first);
    LN* next = iter->next;
    iter->next = map[dest_bin];
    map[dest_bin] = iter;
    rehash_all(next);
}


template<class KEY,class T, int (*thash)(const KEY& a)>
typename HashMap<KEY,T,thash>:: LN** HashMap<KEY,T,thash>::new_hash_table (int bins) {
    LN ** new_h = new LN*[bins];
    for (int i = 0; i < bins; ++i)
        new_h[i] = new LN();
    return new_h;
}



template<class KEY,class T, int (*thash)(const KEY& a)>
void HashMap<KEY,T,thash>::delete_hash_table (LN**& ht, int bins) {
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

template<class KEY,class T, int (*thash)(const KEY& a)>
void HashMap<KEY,T,thash>::Iterator::advance_cursors(){
    while (true){
        if (current.first == ref_map->bins - 1 && current.second->next == nullptr){
            current = Cursor(-1, nullptr);
            return;
        }
        if (current.second->next == nullptr){
            current.second = ref_map->map[++current.first];
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


template<class KEY,class T, int (*thash)(const KEY& a)>
HashMap<KEY,T,thash>::Iterator::Iterator(HashMap<KEY,T,thash>* iterate_over, bool from_begin)
: ref_map(iterate_over), expected_mod_count(ref_map->mod_count) {

    current = (from_begin) ? Cursor(0, ref_map->map[0]) : Cursor(-1, nullptr);
    if (from_begin && current.second->next == nullptr)
        advance_cursors();

}


template<class KEY,class T, int (*thash)(const KEY& a)>
HashMap<KEY,T,thash>::Iterator::~Iterator()
{}


template<class KEY,class T, int (*thash)(const KEY& a)>
auto HashMap<KEY,T,thash>::Iterator::erase() -> Entry {
    if (expected_mod_count != ref_map->mod_count)
        throw ConcurrentModificationError("HashMap::Iterator::erase");
    if (!can_erase)
        throw CannotEraseError("HashMap::Iterator::erase Iterator cursor already erased");
    if (current.first == -1 && current.second == nullptr)
        throw CannotEraseError("HashMap::Iterator::erase Iterator cursor beyond data structure");

    can_erase = false;
    Entry to_return = current.second->value;
    advance_cursors();
    ref_map->erase(to_return.first);
    expected_mod_count = ref_map->mod_count;
    return to_return;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
std::string HashMap<KEY,T,thash>::Iterator::str() const {
    std::ostringstream answer;
    answer << ref_map->str() << "/current=" << current.first  << "/expected_mod_count=" << expected_mod_count << "/can_erase=" << can_erase;
    return answer.str();
}

template<class KEY,class T, int (*thash)(const KEY& a)>
auto  HashMap<KEY,T,thash>::Iterator::operator ++ () -> HashMap<KEY,T,thash>::Iterator& {
    if (expected_mod_count != ref_map->mod_count)
        throw ConcurrentModificationError("HashMap::Iterator::operator ++(int)");

    if (current.first == -1 && current.second == nullptr)
        return *this;

    if (can_erase)
        advance_cursors();
    else
        can_erase = true;

    return *this;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
auto  HashMap<KEY,T,thash>::Iterator::operator ++ (int) -> HashMap<KEY,T,thash>::Iterator {
    if (expected_mod_count != ref_map->mod_count)
        throw ConcurrentModificationError("HashMap::Iterator::operator ++(int)");

    if (current.first == -1 && current.second == nullptr)
        return *this;

    Iterator to_return(*this);
    if (can_erase)
        advance_cursors();
    else
        can_erase = true;

    return to_return;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
bool HashMap<KEY,T,thash>::Iterator::operator == (const HashMap<KEY,T,thash>::Iterator& rhs) const {
    const Iterator* rhsASI = dynamic_cast<const Iterator*>(&rhs);
    if (rhsASI == 0)
        throw IteratorTypeError("HashMap::Iterator::operator ==");
    if (expected_mod_count != ref_map->mod_count)
        throw ConcurrentModificationError("HashMap::Iterator::operator ==");
    if (ref_map != rhsASI->ref_map)
        throw ComparingDifferentIteratorsError("HashMap::Iterator::operator ==");
    return (current == rhsASI->current);
}


template<class KEY,class T, int (*thash)(const KEY& a)>
bool HashMap<KEY,T,thash>::Iterator::operator != (const HashMap<KEY,T,thash>::Iterator& rhs) const {
    return !(*this == rhs);
}


template<class KEY,class T, int (*thash)(const KEY& a)>
pair<KEY,T>& HashMap<KEY,T,thash>::Iterator::operator *() const {
    if (expected_mod_count != ref_map->mod_count)
        throw ConcurrentModificationError("HashMap::Iterator::operator *");
    if (!can_erase || current == Cursor(-1, nullptr)) {
        std::ostringstream where;
        where << " when size = " << ref_map->size();
        throw IteratorPositionIllegal("HashMap::Iterator::operator * Iterator illegal: "+where.str());
    }
    return current.second->value;
}


template<class KEY,class T, int (*thash)(const KEY& a)>
pair<KEY,T>* HashMap<KEY,T,thash>::Iterator::operator ->() const {
    if (expected_mod_count != ref_map->mod_count)
        throw ConcurrentModificationError("HashMap::Iterator::operator ->");
    if (!can_erase || ref_map->used == 0){
        std::ostringstream where;
        where << current.second->value.second << " when size = " << ref_map->size();
        throw IteratorPositionIllegal("HashMap::Iterator::operator -> Iterator illegal: "+where.str());
    }
    return &current.second->value;
}


}

#endif /* HASH_MAP_HPP_ */
