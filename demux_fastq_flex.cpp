#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <atomic>
#include <mutex>

// Mutex to synchronize file writes
std::mutex output_mutex;  
std::atomic<int> total_reads(0);  // Track total reads processed
std::atomic<int> matched_reads(0);  // Track matched reads

// Function to compute reverse complement of a DNA sequence
std::string reverse_complement(const std::string& seq) {
    std::string rev_comp(seq.rbegin(), seq.rend());
    for (char& c : rev_comp) {
        if (c == 'A') c = 'T';
        else if (c == 'T') c = 'A';
        else if (c == 'C') c = 'G';
        else if (c == 'G') c = 'C';
    }
    return rev_comp;
}

// Function to read the tab file and store the indexes with reverse complements
void read_tab_file(const std::string& tab_file, std::unordered_map<std::string, std::pair<std::string, std::string>>& index_dict) {
    std::ifstream infile(tab_file);
    std::string line;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string sample_id, index1, index2;
        std::getline(ss, sample_id, '\t');
        std::getline(ss, index1, '\t');
        std::getline(ss, index2);
        
        // Store the reverse complement of index2 along with index1
        index_dict[sample_id] = {index1, reverse_complement(index2)};
    }
}

// Function to compare two strings with a maximum number of allowed mismatches
bool is_match_with_mismatches(const std::string& seq1, const std::string& seq2, int max_mismatches = 2) {
    if (seq1.size() != seq2.size()) return false;  // sequences must be of the same length
    int mismatches = 0;
    for (size_t i = 0; i < seq1.size(); ++i) {
        if (seq1[i] != seq2[i]) {
            mismatches++;
            if (mismatches > max_mismatches) return false;  // early exit if too many mismatches
        }
    }
    return mismatches <= max_mismatches;
}

// Function to shift sequence by n positions (positive for right, negative for left)
std::string shift_sequence(const std::string& seq, int n) {
    if (n == 0) return seq;
    if (n > 0) {
        return seq.substr(n) + std::string(n, 'N'); // Right shift (add 'N' to the left)
    } else {
        return std::string(-n, 'N') + seq.substr(0, seq.size() + n); // Left shift (add 'N' to the right)
    }
}

// Function to check if two sequences match with frameshift tolerance
bool is_match_with_frameshift_tolerance(const std::string& seq1, const std::string& seq2, int max_mismatches = 1) {
    // Direct match
    if (is_match_with_mismatches(seq1, seq2, max_mismatches)) {
        return true;
    }

    // Shift by 1 to the right and check for mismatches (allowing 1 mismatch)
    if (seq1.size() > 1 && is_match_with_mismatches(shift_sequence(seq1, 1), seq2, max_mismatches)) {
        return true;
    }
    if (seq2.size() > 1 && is_match_with_mismatches(seq1, shift_sequence(seq2, 1), max_mismatches)) {
        return true;
    }

    // Shift by 1 to the left and check for mismatches (allowing 1 mismatch)
    if (seq1.size() > 1 && is_match_with_mismatches(shift_sequence(seq1, -1), seq2, max_mismatches)) {
        return true;
    }
    if (seq2.size() > 1 && is_match_with_mismatches(seq1, shift_sequence(seq2, -1), max_mismatches)) {
        return true;
    }

    // Shift by 2 to the right and check for mismatches (allowing 0 mismatches)
    if (seq1.size() > 2 && is_match_with_mismatches(shift_sequence(seq1, 2), seq2, 0)) {
        return true;
    }
    if (seq2.size() > 2 && is_match_with_mismatches(seq1, shift_sequence(seq2, 2), 0)) {
        return true;
    }

    // Shift by 2 to the left and check for mismatches (allowing 0 mismatches)
    if (seq1.size() > 2 && is_match_with_mismatches(shift_sequence(seq1, -2), seq2, 0)) {
        return true;
    }
    if (seq2.size() > 2 && is_match_with_mismatches(seq1, shift_sequence(seq2, -2), 0)) {
        return true;
    }

    return false;  // No match found
}

// Function to process a pair of reads (record1 and record2) and check for matches
void process_pair(const std::vector<std::string>& record1, const std::vector<std::string>& record2,
                  const std::unordered_map<std::string, std::pair<std::string, std::string>>& index_dict,
                  std::ofstream& unmatched_R1, std::ofstream& unmatched_R2,
                  std::unordered_map<std::string, std::pair<std::ofstream*, std::ofstream*>>& output_files) {

    // Extract indexes from the read description
    size_t idx_start = record1[0].find_last_of(":") + 1;
    std::string read_index = record1[0].substr(idx_start);
    size_t idx_delim = read_index.find("+");
    std::string read_index1 = read_index.substr(0, idx_delim);
    std::string read_index2 = read_index.substr(idx_delim + 1);

    // Iterate over each sample and check for matches
    bool matched = false;
    for (const auto& [sample_id, indexes] : index_dict) {
        const auto& [index1, index2] = indexes;

        // First, try a direct match with allowed mismatches (1 or 2 mismatches)
        if (is_match_with_mismatches(read_index1, index1, 2) && 
            is_match_with_mismatches(read_index2, index2, 2)) {  // Apply 2 mismatch tolerance for both index1 and index2

            matched = true;
            // Write matched records to the corresponding files
            std::lock_guard<std::mutex> lock(output_mutex);
            *output_files.at(sample_id).first << record1[0] << "\n" << record1[1] << "\n" << record1[2] << "\n" << record1[3] << "\n";
            *output_files.at(sample_id).second << record2[0] << "\n" << record2[1] << "\n" << record2[2] << "\n" << record2[3] << "\n";
            matched_reads++;
            break;  // Exit the loop once a match is found
        }

        // If direct match failed, apply frameshift tolerance
        if (is_match_with_frameshift_tolerance(read_index1, index1, 1) && 
            is_match_with_frameshift_tolerance(read_index2, index2, 1)) {  // Apply frameshift comparison for both index1 and index2

            matched = true;
            // Write matched records to the corresponding files
            std::lock_guard<std::mutex> lock(output_mutex);
            *output_files.at(sample_id).first << record1[0] << "\n" << record1[1] << "\n" << record1[2] << "\n" << record1[3] << "\n";
            *output_files.at(sample_id).second << record2[0] << "\n" << record2[1] << "\n" << record2[2] << "\n" << record2[3] << "\n";
            matched_reads++;
            break;  // Exit the loop once a match is found
        }
    }

    // If no match was found, write to unmatched files
    if (!matched) {
        std::lock_guard<std::mutex> lock(output_mutex);
        unmatched_R1 << record1[0] << "\n" << record1[1] << "\n" << record1[2] << "\n" << record1[3] << "\n";
        unmatched_R2 << record2[0] << "\n" << record2[1] << "\n" << record2[2] << "\n" << record2[3] << "\n";
    }
}

// Function to process the FASTQ files in parallel
void process_fastq_files(const std::string& fastq_file_1, const std::string& fastq_file_2,
                         const std::unordered_map<std::string, std::pair<std::string, std::string>>& index_dict,
                         const std::string& output_dir) {

    std::ifstream infile1(fastq_file_1);
    std::ifstream infile2(fastq_file_2);
    std::ofstream unmatched_R1(output_dir + "/unmatched_R1.fastq");
    std::ofstream unmatched_R2(output_dir + "/unmatched_R2.fastq");

    // Create a map for storing output files for each sample
    std::unordered_map<std::string, std::pair<std::ofstream*, std::ofstream*>> output_files;

    // Initialize output files for each sample
    for (const auto& [sample_id, indexes] : index_dict) {
        output_files[sample_id] = {new std::ofstream(output_dir + "/" + sample_id + "_R1.fastq"),
                                   new std::ofstream(output_dir + "/" + sample_id + "_R2.fastq")};
    }

    std::string line1, line2;
    while (std::getline(infile1, line1) && std::getline(infile2, line2)) {
        std::vector<std::string> record1{line1, "", "", ""};
        std::vector<std::string> record2{line2, "", "", ""};

        // Read next 3 lines for each read
        for (int i = 1; i < 4; ++i) {
            std::getline(infile1, record1[i]);
            std::getline(infile2, record2[i]);
        }

        total_reads++;
        process_pair(record1, record2, index_dict, unmatched_R1, unmatched_R2, output_files);
    }

    // Close all output files
    for (auto& [sample_id, files] : output_files) {
        files.first->close();
        files.second->close();
        delete files.first;
        delete files.second;
    }

    unmatched_R1.close();
    unmatched_R2.close();
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <tab_file> <fastq_file_1> <fastq_file_2> <output_dir>\n";
        return 1;
    }

    std::string tab_file = argv[1];
    std::string fastq_file_1 = argv[2];
    std::string fastq_file_2 = argv[3];
    std::string output_dir = argv[4];

    // Load index information from the tab file
    std::unordered_map<std::string, std::pair<std::string, std::string>> index_dict;
    read_tab_file(tab_file, index_dict);

    // Process FASTQ files
    process_fastq_files(fastq_file_1, fastq_file_2, index_dict, output_dir);

    std::cout << "Total reads processed: " << total_reads.load() << std::endl;
    std::cout << "Matched reads: " << matched_reads.load() << std::endl;
    std::cout << "Unmatched reads written to unmatched_R1.fastq and unmatched_R2.fastq" << std::endl;

    return 0;
}
