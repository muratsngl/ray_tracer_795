#!/bin/bash

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Counters
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Arrays to store results
PASSED_FILES=()
FAILED_FILES=()

echo -e "${BLUE}==========================================${NC}"
echo -e "${BLUE}  Ray Tracer Scene Parser Test Suite${NC}"
echo -e "${BLUE}==========================================${NC}"
echo ""

# Make sure the raytracer is built
echo -e "${YELLOW}Building raytracer...${NC}"
if make; then
    echo -e "${GREEN}✓ Build successful${NC}"
else
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
fi
echo ""

# Function to test a single JSON file
test_scene() {
    local json_file="$1"
    local relative_path="${json_file#*/inputs/}"
    
    echo -n "Testing $relative_path... "
    
    # Run the raytracer and capture output
    output=$(timeout 30 ./raytracer "$json_file" 2>&1)
    exit_code=$?
    
    # Check if the test passed
    if [ $exit_code -eq 0 ]; then
        echo -e "${GREEN}✓ PASSED${NC}"
        PASSED_FILES+=("$relative_path")
        ((PASSED_TESTS++))
        
        # Extract and display key parsing info
        vertex_count=$(echo "$output" | grep "Vertices:" | awk '{print $2}')
        material_count=$(echo "$output" | grep "Materials:" | awk '{print $2}')
        camera_count=$(echo "$output" | grep "Cameras:" | awk '{print $2}')
        
        if [[ -n "$vertex_count" ]]; then
            echo -e "    ${BLUE}→${NC} Vertices: $vertex_count, Materials: $material_count, Cameras: $camera_count"
        fi
    elif [ $exit_code -eq 124 ]; then
        echo -e "${RED}✗ TIMEOUT${NC}"
        FAILED_FILES+=("$relative_path (TIMEOUT)")
        ((FAILED_TESTS++))
    else
        echo -e "${RED}✗ FAILED${NC}"
        FAILED_FILES+=("$relative_path")
        ((FAILED_TESTS++))
        
        # Show error details
        echo -e "    ${RED}Error:${NC} $output" | head -n 3
    fi
    
    ((TOTAL_TESTS++))
    echo ""
}

# Find and test all JSON files
echo -e "${YELLOW}Discovering JSON scene files...${NC}"

# Test main directory JSON files
for json_file in inputs/*.json; do
    if [ -f "$json_file" ]; then
        test_scene "$json_file"
    fi
done

# Test subdirectory JSON files
for json_file in inputs/*/*.json; do
    if [ -f "$json_file" ]; then
        test_scene "$json_file"
    fi
done

# Test nested subdirectory JSON files (for raven folder)
for json_file in inputs/*/*/*.json; do
    if [ -f "$json_file" ]; then
        test_scene "$json_file"
    fi
done

# Print summary
echo -e "${BLUE}==========================================${NC}"
echo -e "${BLUE}           TEST SUMMARY${NC}"
echo -e "${BLUE}==========================================${NC}"
echo -e "Total tests run: ${BLUE}$TOTAL_TESTS${NC}"
echo -e "Tests passed:    ${GREEN}$PASSED_TESTS${NC}"
echo -e "Tests failed:    ${RED}$FAILED_TESTS${NC}"

if [ $PASSED_TESTS -eq $TOTAL_TESTS ]; then
    echo -e "\n${GREEN}🎉 ALL TESTS PASSED! 🎉${NC}"
else
    echo -e "\nSuccess rate: $(( PASSED_TESTS * 100 / TOTAL_TESTS ))%"
fi

# List passed files
if [ ${#PASSED_FILES[@]} -gt 0 ]; then
    echo -e "\n${GREEN}✓ Passed files:${NC}"
    for file in "${PASSED_FILES[@]}"; do
        echo -e "  ${GREEN}→${NC} $file"
    done
fi

# List failed files
if [ ${#FAILED_FILES[@]} -gt 0 ]; then
    echo -e "\n${RED}✗ Failed files:${NC}"
    for file in "${FAILED_FILES[@]}"; do
        echo -e "  ${RED}→${NC} $file"
    done
    echo ""
    echo -e "${YELLOW}Tip: Run individual files manually for detailed error info:${NC}"
    echo -e "  ${BLUE}./raytracer inputs/filename.json${NC}"
fi

echo ""
echo -e "${BLUE}==========================================${NC}"

# Exit with appropriate code
if [ $FAILED_TESTS -eq 0 ]; then
    exit 0
else
    exit 1
fi