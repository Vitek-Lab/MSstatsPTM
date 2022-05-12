#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
CharacterVector extract_protein_name(CharacterVector ptm_list,
                                     CharacterVector protein_list){
  
  int ptm_len = ptm_list.length();
  int protein_len = protein_list.length();
  CharacterVector identified_proteins;
  
  for (int i=0;i < ptm_len; i++){
    
    // No match is NA
    String protein_match = "NA";
    
    for (int j=0;j < protein_len; j++){
      
      // Initialize temp variables
      String tracker = ptm_list[i];
      String temp_protein = protein_list[j];
      String temp_ptm = ptm_list[i];
      temp_ptm.replace_first(temp_protein, "");
      
      // Check for match
      if (temp_ptm != tracker){
        // Take first match, sorted by string size before Rcpp function
        protein_match = temp_protein;
        break;
      };
    };
    identified_proteins.push_back(protein_match);
  };
  return identified_proteins;
}
