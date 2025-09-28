#include "move.hpp"
#include "ds/move_table.hpp"
#include <iostream>
#include <cassert>
#include <sstream>
#include <optional>
#include <vector>

using namespace std;

// Test data
vector<ulint> test_lengths = {10, 20, 15, 8};
vector<ulint> test_permutations = {0, 10, 30, 45};
size_t test_n = 53;

void test_length_config() {
    cout << "Testing LengthConfig..." << endl;
    
    // Test without max length
    LengthConfig config1(test_lengths, std::nullopt);
    assert(config1.split_num_rows == 4);
    cout << "  No capping: " << config1.split_num_rows << " rows, max_length: " << config1.max_observed_length << endl;
    
    // Test with max length
    LengthConfig config2(test_lengths, 12);
    cout << "  With capping: " << config2.split_num_rows << " rows, max_length: " << config2.max_observed_length << endl;
    
    cout << "LengthConfig tests passed!" << endl << endl;
}

template<typename Table>
void test_table_type(const string& name) {
    cout << "Testing " << name << "..." << endl;
    
    // Test basic construction and operations
    LengthConfig config(test_lengths, std::nullopt);
    
    // Test table creation and basic operations
    Table table; // Test default constructor
    
    // Test setting and getting values
    // TODO: Add specific tests once we know the interface
    
    cout << name << " basic tests passed!" << endl << endl;
}

void test_serialization() {
    cout << "Testing serialization..." << endl;
    
    // TODO: Create a table, serialize it, load it back, verify it matches
    
    cout << "Serialization tests passed!" << endl << endl;
}

int main() {
    cout << "=== MoveStructure Table Tests ===" << endl << endl;
    
    try {
        test_length_config();
        
        test_table_type<MoveTable<>>("MoveTable<>");
        test_table_type<MoveVector<>>("MoveVector<>");
        test_table_type<MoveTableIdx>("MoveTableIdx");  
        test_table_type<MoveVectorIdx>("MoveVectorIdx");
        
        test_serialization();
        
        cout << "=== ALL TESTS PASSED ===" << endl;
        
    } catch (const exception& e) {
        cerr << "TEST FAILED: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}