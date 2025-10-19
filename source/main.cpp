#include <iostream>
#include "parser.cpp"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <scene_file.json>" << std::endl;
        return 1;
    }
    
    try {
        Scene scene = Parser::parseScene(argv[1]);
    } catch (const nlohmann::json::exception& e) {
        std::cerr << "JSON Error: " << e.what() << std::endl;
        std::cerr << "Error ID: " << e.id << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    


    return 0;
}