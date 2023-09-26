#include "shortest_certificate.h"
#include "frechet_light.h"
#include "freespace_light_vis.h"

#include <map>
#include <fstream>

void printUsage() {
	std::cout <<
		"Usage: ./shortest_certificate_bench <characters | sigspatial | geolife> <yes | yes_alt | no>\n"
		"\n";
}

int main(int argc, char* argv[]) {

    if(argc != 3) {
        printUsage();
        return 1;
    }

    std::string file_path;
    std::string file_dest;
    bool cert_type = (std::strcmp(argv[2], "no") == 0) ? 0 : 1;

    if(std::strcmp(argv[1], "characters") == 0) {
        std::cout << "Executing characters " << argv[2] << "..." << std::endl;
        if(cert_type) file_path = "../test_data/decider_benchmark_queries/characters_query_decider_1_0_plus.txt";
        else file_path = "bench_data/query_characters_no.txt";
        file_dest = "../test_data/benchmark/characters/data/";
    }
    else if(std::strcmp(argv[1], "sigspatial") == 0) {
        std::cout << "Executing sigspatial " << argv[2] << "..." << std::endl;
        if(cert_type) file_path = "../test_data/decider_benchmark_queries/sigspatial_query_decider_2_-5_plus.txt";
        else file_path = "bench_data/query_sigspatial_no.txt";
        file_dest = "../test_data/benchmark/sigspatial/";
    }
    else if(std::strcmp(argv[1], "geolife") == 0) {
        std::cout << "Executing geolife " << argv[2] << "..." << std::endl;
        if(cert_type) file_path = "../test_data/decider_benchmark_queries/geolife_query_decider_2_-4_plus.txt";
        else file_path = "bench_data/query_geolife_no.txt";
        file_dest = "../test_data/benchmark/geolife/Data/";
    }
    else {
        std::cout << "Wrong argument..." << std::endl;
        return 1;
    } 

    std::string argv1 = argv[1];
    std::string argv2 = argv[2];

    std::map<uint32_t, uint32_t> results;
    std::map<uint32_t, uint32_t> results_compare;

    std::ifstream file(file_path);
    std::ofstream out("data/" + argv1 + "_" + argv2 + ".txt");
    std::string c1, c2, delta;
    uint32_t counter = 0;

    while(file >> c1 >> c2 >> delta) {
        ShortestCertificate s(file_dest + c1, file_dest + c2, delta);
        FrechetLight frechet;
        Certificate result;
	
        //Bound for cubic runtime
        if(s.curve1.size() < 1000 && s.curve2.size() < 1000 && s.curve1.size()) {
            counter++;
            std::cout << "[" << s.curve1.size()<< "," << s.curve2.size() <<"]" << std::endl;

            std::cout << "FRECHET_LIGHT:\n";
            auto timer_start = std::chrono::high_resolution_clock::now();
            bool output = frechet.lessThan(std::stod(delta), s.curve1, s.curve2);
	        Certificate& c = frechet.computeCertificate();
            c.dump_certificate();
            auto time = std::chrono::high_resolution_clock::now() - timer_start;
            auto time_in_s = std::chrono::duration_cast<std::chrono::milliseconds>(time);           
            std::cout << "["<< (counter) << "]" << "[" << time_in_s.count() << "ms]\n\n";

            std::cout << "SHORTEST_CERTIFICATE:\n";
            timer_start = std::chrono::high_resolution_clock::now();
            if(std::strcmp(argv[2], "yes_alt") == 0) result = s.yes_certificate(s.curve1, s.curve2, s.matrix, s.delta);
            else result = cert_type ? s.yes_certificate(s.curve1, s.curve2, s.matrix, s.delta) : s.no_certificate(s.curve1, s.curve2, s.matrix, s.delta);         
            time = std::chrono::high_resolution_clock::now() - timer_start;
            time_in_s = std::chrono::duration_cast<std::chrono::milliseconds>(time);           
            std::cout << "["<< (counter) << "]" << "[" << time_in_s.count() << "ms]\n\n";

            /*if(result.size() == 1 || result.size() == 0) {
                continue;
            }else counter++;*/
            
            //result.check();
            c.check();

            std::cout << c1 << " " << c2 << " " << delta << std::endl;
            
            auto it = results.find(result.size());
            if(it != results.end()) {
                it->second += 1;
            }
            else results.insert(std::make_pair(result.size(), 1));

            it = results_compare.find(c.size());
            if(it != results_compare.end()) {
                it->second += 1;
            }
            else results_compare.insert(std::make_pair(c.size(), 1));
                      
        }

        /*if(counter == 0) {
            //frechet.setCertificate(result);
            FreespaceLightVis vis(frechet);
            vis.exportFreespaceToSvg("data/plot_compare_yes_fl_adap.svg");
            frechet.setCertificate(result);
            FreespaceLightVis vis2(frechet);
            vis2.exportFreespaceToSvg("data/plot_diagonal_yes.svg");
        }*/

        if(counter == 1000) break;
    }
    file.close();
    

    for(auto elem : results) {
        std::cout << elem.first << " " << elem.second << std::endl;
        out << elem.first << " " << elem.second << std::endl;
    }
    out << "-" << std::endl;
    std::cout << "\n\n";
    for(auto elem : results_compare) {
        std::cout << elem.first << " " << elem.second << std::endl;
        out << elem.first << " " << elem.second << std::endl;
    }
    out << "-" << std::endl;
    out.close();
}
