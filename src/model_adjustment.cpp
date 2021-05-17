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


// int rstring_length( SEXP sv ){
//   SEXP s = STRING_ELT(sv, 0) ;
//   return strlen( CHAR(s) ) ;
// }
// 
// 
// 
// DataFrame cpp_rbind(List df_list){
// 
//   int n_lists = df_list.length();
//   int expand_rows = 0;
//   for (int i=0;i < n_lists;i++){
//     expand_rows += df_list[i].length();
//   };
// 
//   CharacterVector Protein = CharacterVector(expand_rows);
//   CharacterVector Site = CharacterVector(Site);
//   CharacterVector Label = CharacterVector(Label);
//   NumericVector log2FC = NumericVector(log2FC);
//   NumericVector SE = NumericVector(SE);
//   NumericVector Tvalue = NumericVector(Tvalue);
//   NumericVector DF = NumericVector(DF);
//   NumericVector Pvalue = NumericVector(Pvalue);
//   NumericVector Adj_pvalue = NumericVector(Adj_pvalue);
// 
//   int start = 0;
//   for (int i=0;i < n_lists;i++){
//     DataFrame temp_df = df_list[i];
//     int temp_length = temp_df.length()
//     Protein[start:temp_length] = temp_df["Protein"];
//     Site[start:temp_length] = temp_df["Site"];
//     Label[start:temp_length] = temp_df["Label"];
//     log2FC[start:temp_length] = temp_df["log2FC"];
//     SE[start:temp_length] = temp_df["SE"];
//     Tvalue[start:temp_length] = temp_df["Tvalue"];
//     DF[start:temp_length] = temp_df["DF"];
//     Pvalue[start:temp_length] = temp_df["pvalue"];
//     Adj_pvalue[start:temp_length] = temp_df["adj.pvalue"];
// 
//     start += temp_length;
//     }
// 
//   return DataFrame::create(Named("Protein") = Protein,
//                            Named("Site") = Site,
//                            Named("Label") = Label,
//                            Named("log2FC") = log2FC,
//                            Named("SE") = SE,
//                            Named("Tvalue") = Tvalue,
//                            Named("DF") = DF,
//                            Named("Pvalue") = Pvalue,
//                            Named("Adj_pvalue") = Adj_pvalue);
// 
// }

// DataFrame apply_ptm_adjustment(DataFrame ptm_model,
//                                DataFrame protein_model,
//                                CharacterVector comparison_vector) {
//
//   int n_comparison = comparison_vector.length();
//
//   for (int i=0;i<n_comparison;i++){
//
//     // Filter for label
//     DataFrame sub_ptm = filter_dataframe(ptm_model, comparison_vector[i]);
//     DataFrame sub_protein = filter_dataframe(protein_model,
//                                                 comparison_vector[i]);
//
//   }
//   return sub_protein
// }
//
// DataFrame filter_dataframe(DataFrame df,
//                            CharacterVector value) {
//   // Filter for label
//   StringVector sub = df["Label"];
//   LogicalVector ind(sub.size());
//   for (int i=0; i<sub.size();i++){
//     ind[i] = (sub(i) == value[0]);
//   }
//
//   CharacterVector Protein = df["Protein"];
//   CharacterVector Site = df["Site"];
//   CharacterVector Label = df["Label"];
//   NumericVector log2FC = df["log2FC"];
//   NumericVector SE = df["SE"];
//   NumericVector Tvalue = df["Tvalue"];
//   NumericVector DF = df["DF"];
//   NumericVector pvalue = df["pvalue"];
//
//   return DataFrame::create(Named("Protein") = Protein[ind],
//                            Named("Site") = Site[ind],
//                            Named("Label") = sub[ind],
//                            Named("log2FC") = log2FC[ind],
//                            Named("SE") = SE[ind],
//                            Named("Tvalue") = Tvalue[ind],
//                            Named("DF") = DF[ind],
//                            Named("pvalue") = pvalue[ind]
//                           );
// }


