# Context Free Grammar functions
# Moritz Feigl, 2019
#


rule <- function(rule_string) {
  # Helper function for grammar definition
  
  # split by comma
  rule <- unlist(strsplit(rule_string, split = ", "))
  # find parantheses and fix possible wrong splits
  ind_par_open <- grep(pattern = "(", rule, fixed = TRUE)
  if(length(ind_par_open) != 0){
    ind_par_close <- grep(pattern = ")", rule, fixed = TRUE)
    if(length(ind_par_open) != length(ind_par_close)){
      stop("Rule formulation error. Single parenthesis found.")
    }
    del <- integer()
    for(i in seq_along(ind_par_open)){
      if(ind_par_open[i] != ind_par_close[i]){
        rule[ind_par_open[i]] <- paste0(rule[ind_par_open[1]], rule[ind_par_close[i]])
        del <- c(del, ind_par_close[i])
      }
    }
    # delete wrong splits
    if(length(del != 0)) rule <- rule[-del]
  }
  
  return(rule)
}

Grammar <- setClass(
  # Set the name for the class
  "Grammar",
  
  # Define the slots
  slots = c(
    input = "list",
    grammar = "list"
  )
)

grammar <- function(...) {
  ### Create a Grammar object from a list of rules ###
  
  def <- list(...)
  # rename non-terminals
  names(def) <- paste0("<", names(def), ">")
  # apply rule function
  def_gram <- lapply(def, rule)
  # create class Grammar
  def_gram <- Grammar(input = def,
                      grammar = def_gram)
  return(def_gram)
}


setMethod("show", "Grammar", function(object) {
  ### print method for Grammar object ###
  print(object@input)
  
})


grammar.unroller <- function(Grammar, max.depth){
  ### unrolls a recursive grammar ###
  
  # go through all non-terminals and roll the recursive non-terminals out to max depth
  unrolled_grammar <- lapply(1:length(Grammar), nt.unroller, Grammar = Grammar, max.depth = max.depth)
  unrolled_grammar_unique <- vector(mode = "list")
  # make each recursive non-terminal a named element of the list
  for(k in 1:length(unrolled_grammar)){
    gram_part <- unrolled_grammar[[k]]
    # if grammar part is larger than one row -> split
    if(!is.null(nrow(gram_part))){
      gram_part_list <- vector(mode = "list", length = nrow(gram_part))
      names(gram_part_list) <- rownames(gram_part)
      for (i in 1:nrow(gram_part)) gram_part_list[[i]] <- gram_part[i, ]
      unrolled_grammar_unique <- c(unrolled_grammar_unique, gram_part_list)
    } else {
      unrolled_grammar_unique <- c(unrolled_grammar_unique, unrolled_grammar[k])
    }
  }
  return(unrolled_grammar_unique)
}

nt.unroller <- function(Grammar, nt_id, max.depth){
  ### unrolls terminal to max depth ### 
  
  # get NT name and statements
  nt <- names(Grammar)[nt_id]
  nts <- Grammar[[nt_id]]
  # is recursive?
  if(sum(grepl(nt, nts)) > 0) {
    rnts <- data.frame()
    # Recursive unrolling
    for(i in 1:(max.depth + 1)) rnts <- rbind(rnts, nts, stringsAsFactors = FALSE)
    # Make changed non-terminal at max depth:
    rnts <- rbind(rnts, gsub("<.*?>", "NA", rnts[1, ]), stringsAsFactors = FALSE)
    # check last row and change to usefull arguments
    last_row <- rnts[nrow(rnts), ]
    last_row <- as.character(last_row)
    for(i in seq_along(last_row)){
      tmp <- unlist(strsplit(last_row[i], "+", fixed = TRUE))
      tmp <- unlist(strsplit(tmp, "-", fixed = TRUE))
      if(length(tmp) < 2) {
        tmp <- "NA"
        last_row[i] <- tmp
      }
    }
    rnts[nrow(rnts), ] <- last_row
    # col and rownames
    colnames(rnts) <- NULL
    rownames(rnts) <- c(nt, paste0(substr(nt, 1, (nchar(nt)-1) ), 1:(max.depth + 1), ">"))
    # change names corresponding to new non-terminals
    for(i in 1:(nrow(rnts) - 1)){
      rnts[i, ] <- gsub(rownames(rnts)[1], rownames(rnts)[i+1], rnts[i, ])
    }
    return(rnts)
  } else {
    # make a data.frame to define rownames
    nts <- data.frame(t(nts), stringsAsFactors = FALSE)
    colnames(nts) <- NULL
    rownames(nts) <- nt
    return(nts)
  }
}

na.remover <- function(grammar_function){
  ### removes NAs from unrolled grammar ###
  
  # remove all spaces
  gf <- gsub(" ", "", grammar_function)
  # split everything
  gf_split <- unlist(strsplit(gf, "+"))
  # reassable NA strings
  for(i in seq_along(gf_split[-1])){
    if(gf_split[i] == "N" & gf_split[i+1] == "A") gf_split[i] <- "NA"
  }
  # if there are no more NAs return
  if(sum(gf_split == "NA") == 0) return(paste0(gf_split, collapse = ""))
  # remove NAs
  gf_split <- gf_split[-(which(gf_split == "NA") + 1)]
  # if there are some NAs just next to each other -> combine into one
  nas <- which(gf_split == "NA")
  while(length(nas) > 1) {
    if(nas[2] == (nas[1] + 1)) {
      gf_split <- gf_split[-nas[1]]
      nas <- nas[-1]
    } else {
      nas <- nas[-1]
    }
  }
  infinite_detector <- 0
  while(sum(gf_split == "NA") != 0){
    # catch infinite recursion
    infinite_detector <- infinite_detector + 1
    if (infinite_detector > length(gf_split)){
      stop("Infinite recursion detected. Check if there are any typos in your Grammar. There might be a terminal adjacent to another without a +, -, *, /")
    }
    # if there is a NA function of the form NA(...) replace it with NA
    na_ind <- which(gf_split == "NA")
    delete <- matrix(nrow = 1, ncol = 2)
    for (i in na_ind){
      if(i != length(gf_split)){
        # the last NA is left out -> can't be a function
        if(gf_split[i + 1] == "("){
          # find the following ")"
          fclose_ind <- which(gf_split == ")")
          fclose_ind <- fclose_ind[fclose_ind > i][1]
          # position i stays NA and rest will afterwards be deleted
          delete <- rbind(delete, c(i, fclose_ind))
        }
      }
    }
    if(nrow(delete) > 1){
      # delete previously chosen entries
      for(i in nrow(delete):2){
        gf_split <- gf_split[-((delete[i, 1]+1):delete[i, 2])]
      }
    }
    # added or substracted -> NAs are turned to 0
    na_ind <- which(gf_split == "NA")
    for (i in na_ind){
      if(i == 1){
        if(gf_split[i + 1] %in% c("+", "-")) gf_split[i] <- 0
      } else if(i == length(gf_split)){
        if(gf_split[i -1] %in% c("+", "-")) gf_split[i] <- 0
      } else {
        if(gf_split[i - 1] %in% c("+", "-") & gf_split[i + 1] %in% c("+", "-") |
           gf_split[i - 1] %in% c("+", "-") & gf_split[i + 1] == ")" |
           gf_split[i - 1] == "(" & gf_split[i + 1] %in% c("+", "-")
        ){
          gf_split[i] <- 0
        }
      }
    }
    # if there is a * or / -> NAs are turned into 1
    if(sum(gf_split == "NA") == 0) return(paste0(gf_split, collapse = ""))
    na_ind <- which(gf_split == "NA")
    for (i in na_ind){
      if(i == 1){
        if(gf_split[i + 1] %in% c("*", "/")) gf_split[i] <- 1
      } else if(i == length(gf_split)){
        if(gf_split[i -1] %in% c("*", "/")) gf_split[i] <- 1
      } else {
        if(gf_split[i - 1] %in% c("*", "/") | gf_split[i + 1] %in% c("*", "/")){
          gf_split[i] <- 1
        }
      }
    }
    # if NA is in a function -> turn into 1 as well
    if(sum(gf_split == "NA") == 0) return(paste0(gf_split, collapse = ""))
    na_ind <- which(gf_split == "NA")
    for (i in na_ind){
      if(i != 1){
        if(gf_split[i - 1] == "(" & gf_split[i + 1] == ")") gf_split[i] <- 1
      }
    }
  }
  return(paste0(gf_split, collapse = ""))
}

search.table <- function(Grammar, max.depth = NULL){
  ### creates a vector representation of a given Grammar ###
  
  # subset s4 object
  if(!is.list(Grammar)){
    Grammar <- Grammar@grammar
  }
  # remove space from grammar
  Grammar <- lapply(Grammar, function(x) gsub(" ", "", x))
  # Recursive or not?
  nonterm <- names(Grammar)
  rec <- vector(mode ="list", length = length(Grammar))
  for(i in 1:length(Grammar)){
    rec[[i]] <- grep(nonterm[i], Grammar[[i]])
  }
  # Which are the recursive non-terminals
  r_nt <- names(Grammar)[sapply(rec, function(x) length(x) > 0)]
  # Do the recursive non-terminals include other recursive non-terminals?
  # loop over all recursive non-terminals to investigate if they have other rec. nts
  r_nt_pointer <- vector(mode = "list", length = length(r_nt))
  names(r_nt_pointer) <- r_nt
  for(l in seq_along(r_nt)){
    # loop over all recursive non-terminals that are not investigated
    for(n in seq_along(r_nt[-l])) {
      is_rec.nt_inrec.nt <- grepl(r_nt[-l][n], Grammar[r_nt[l]])
      if(is_rec.nt_inrec.nt) r_nt_pointer[[l]] <- c(r_nt_pointer[[l]], r_nt[-l][n])
    }
  }
  # if double recursive elements are present
  if(length(r_nt_pointer) > 0){
    # remove non-terminals that have no other recursive non-terminals included in their statement
    r_nt_pointer <- r_nt_pointer[!unlist(lapply(r_nt_pointer, is.null))]
  }
  # find out if there is one recursive nt that points to more than one recursive nt that points back to it
  r_nt_names <- names(r_nt_pointer)
  finding_deep_recursions <- r_nt_pointer
  for(i in seq_along(r_nt_names)){
    for(j in 1:length(r_nt_pointer[[i]]))
      finding_deep_recursions[[i]][j] <- sum(grepl(r_nt_names[i], r_nt_pointer[[r_nt_pointer[[i]][j]]])) > 0
  }
  # if any recursive nt hase more than one recursive nt pointing at him -> error
  for(i in seq_along(r_nt_names)){
    if(sum(as.logical(finding_deep_recursions[[i]])) > 1){
      stop("Recursion width is too deep. At most, only two recursive non-terminals are allowed to refer to each other.")
    }
  }
  # now r_nt_pointer is a named list that includes all recrusive non-terminals that
  # point to another recursive non-terminal. The names are the non-terminals, and the
  # respective list elements are the other rec. non-terminals they point to.
  # if recursive, apply unroll function
  if (length(r_nt) > 0){
    # Recursive search table ->
    if(is.null(max.depth)) stop("Grammar is recursive. Please supply a maximum depth length.")
    # unroll Grammar to include recursive terms, depending on the max.depth
    urGrammar <- grammar.unroller(Grammar, max.depth = max.depth)
  } else {
    urGrammar <- Grammar
  }
  # double recursion:
  # in case there is a double recursion, the higher non-terminal will be set to the value at max.depth
  # define Grammar hierachy
  gHierachy <- data.frame(nt = names(Grammar), hierachy = 1:length(Grammar))
  # Do only if double recursions exist
  r_nt_double_rec <- vector(mode = "list", length(r_nt_pointer))
  if(length(r_nt_pointer) > 0){
    names(r_nt_double_rec) <- names(r_nt_pointer)
    for (i in 1:length(r_nt_pointer)){
      NTi <- r_nt_pointer[[i]]
      NTi_name <- names(r_nt_pointer)[i]
      # is one of the pointed recursive elements lower in the hierachy (later in Grammar)?
      hierachy_of_nt <- gHierachy[gHierachy$nt == NTi_name, "hierachy"]
      hierachy_of_pointer <- gHierachy[gHierachy$nt %in% NTi, "hierachy"]
      r_nt_double_rec[[i]] <- NTi[hierachy_of_nt > hierachy_of_pointer]
    }
    # remove non-terminals that have no double recursion
    r_nt_double_rec <- r_nt_double_rec[unlist(lapply(r_nt_double_rec, function(x) length(x) > 0))]
  }
  # if there are still double recursive elements
  if(length(r_nt_double_rec) > 0){
    # change double recursion links: all that point upwards to a recursion are set to max.depth
    number_of_changes <- numeric()
    for(i in 1:length(r_nt_double_rec)){
      nt_to_change <- names(r_nt_double_rec)[i]
      nt_to_change_into <- r_nt_double_rec[[i]]
      nts_to_change <- names(urGrammar)[grep(nt_to_change, gsub("[1-9]", "", names(urGrammar)))]
      for(k in seq_along(nts_to_change)){
        new_nt_name <- paste0(substr(nt_to_change_into, 1, nchar(nt_to_change_into) - 1), max.depth + 1, ">")
        urGrammar[[nts_to_change[k]]] <- gsub(nt_to_change_into, new_nt_name, urGrammar[[nts_to_change[k]]])
      }
      number_of_changes <- c(number_of_changes, length(nts_to_change))
    }
    # produce list with added number of max.depth objects due to double recursion
    for(i in 1:length(r_nt_double_rec)) {
      r_nt_double_rec[[i]] <- cbind(r_nt_double_rec[[i]], number_of_changes[i])
      r_nt_double_rec[[i]] <- cbind(r_nt_double_rec[[i]],
                                    paste0(substr(r_nt_double_rec[[i]][1], 1,
                                                  nchar(r_nt_double_rec[[i]][1]) - 1),
                                           max.depth + 1,
                                           ">"))
    }
  }
  # Double pointing without recursion
  #  Find out if there are any nts that point at each other and create a infinite recursion.
  #  Don't take in those that are already dealt with in the double recursive part.
  nt_pointer <- vector(mode = "list", length = length(Grammar))
  names(nt_pointer) <- nonterm
  if(length(r_nt_pointer) > 0) {
    nt_pointer <- nt_pointer[!(names(nt_pointer) %in% r_nt_pointer)]
  }
  pointer <- names(nt_pointer)
  # remove nts that have no pointers inside
  ind_empty_point <- logical(length = length(pointer))
  for(l in seq_along(nt_pointer)){
    ind_empty_point[l] <- grepl("<*?>", Grammar[pointer[l]])
  }
  # remove those
  nt_pointer <- nt_pointer[ind_empty_point]
  pointer <- pointer[ind_empty_point]
  for(l in seq_along(nt_pointer)){
    # loop over all non-terminals that are not investigated
    for(n in seq_along(pointer[-l])) {
      points_to <- grepl(pointer[-l][n], Grammar[pointer[l]])
      points_back <- grepl(pointer[l], Grammar[pointer[-l][n]])
      
      if(sum(points_to, points_back) == 2) nt_pointer[[l]] <- c(nt_pointer[[l]], pointer[-l][n])
    }
  }
  # remove empty nt_pointer
  if(length(nt_pointer) > 0) nt_pointer <- nt_pointer[!unlist(lapply(nt_pointer, is.null))]
  if(length(nt_pointer) > 0){
    # double pointer processing needs a max.depth
    if(is.null(max.depth)) stop("Grammar is recursive. Please supply a maximum depth length.")
  }
  # in case there is more than a double pointing interaction -> stop function
  nt_pointers_larger1 <- unlist(lapply(nt_pointer, function(x) length(x) > 1))
  if(sum(nt_pointers_larger1) > 0){
    stop("Recursion width is too deep. At most, only two non-terminals are allowed to refer to each other.")
  }
  # in case there is a double pointer, the higher non-terminal will be set to the value at max.depth
  # Do only if double pointer exist
  nt_double_pointer <- vector(mode = "list", length(nt_pointer))
  if(length(nt_pointer) > 0){
    names(nt_double_pointer) <- names(nt_pointer)
    for (i in seq_along(nt_pointer)){
      NTi <- nt_pointer[[i]]
      NTi_name <- names(nt_pointer)[i]
      # is one of the pointed recursive elements lower in the hierachy (later in Grammar)?
      hierachy_of_nt <- gHierachy[gHierachy$nt == NTi_name, "hierachy"]
      hierachy_of_pointer <- gHierachy[gHierachy$nt %in% NTi, "hierachy"]
      nt_double_pointer[[i]] <- NTi[hierachy_of_nt > hierachy_of_pointer]
    }
    # remove non-terminals that are higher in hierachy that the ones they point to
    nt_double_pointer <- nt_double_pointer[unlist(lapply(nt_double_pointer, function(x) length(x) > 0))]
    
    # change double links: all that point upwards to a recursion are set to depth.id
    if(length(nt_double_pointer) > 0){
      number_of_changes <- numeric()
      for(i in 1:length(nt_double_pointer)){
        nt_to_change <- names(nt_double_pointer)[i]
        nt_to_change_into <- nt_double_pointer[[i]]
        nts_to_change <- names(urGrammar)[grep(nt_to_change, gsub("[1-9]", "", names(urGrammar)))]
        # get highest number of current Grammar for nts_to_change_into
        search_nt_to_change_into <- substr(nt_to_change_into, 1, nchar(nt_to_change_into) - 1)
        max.change <- names(urGrammar)[grep(search_nt_to_change_into, names(urGrammar))]
        max.change <- gsub(search_nt_to_change_into, "", max.change)
        max.change <- gsub(">", "", max.change)
        if(sum(max.change %in% "") != length(max.change)){
          depth.id <- max(as.numeric(max.change), na.rm = TRUE)
        } else depth.id <- 1
        
        for(k in seq_along(nts_to_change)){
          new_nt_name <- paste0(substr(nt_to_change_into, 1, nchar(nt_to_change_into) - 1), depth.id, ">")
          urGrammar[[nts_to_change[k]]] <- gsub(nt_to_change_into, new_nt_name, urGrammar[[nts_to_change[k]]])
        }
        number_of_changes <- c(number_of_changes, length(nts_to_change))
      }
      # produce list with added number of max.depth objects due to double recursion
      for(i in 1:length(nt_double_pointer)) {
        nt_double_pointer[[i]] <- cbind(nt_double_pointer[[i]], number_of_changes[i])
        nt_double_pointer[[i]] <- cbind(nt_double_pointer[[i]],
                                        paste0(substr(nt_double_pointer[[i]][1], 1,
                                                      nchar(nt_double_pointer[[i]][1]) - 1),
                                               depth.id,
                                               ">"))
      }
      # Add new terms max.depth terms to Grammar
      nt_double_pointer_df <- do.call(rbind, nt_double_pointer)
      for(i in 1:nrow(nt_double_pointer_df)){
        new.nt <- vector(mode = "list", length = 1)
        names(new.nt) <- nt_double_pointer_df[i, 3]
        new.nt[[1]] <- urGrammar[[nt_double_pointer_df[i, 1]]]
        # remove all nts
        new.nt[[1]] <- gsub(names(nt_double_pointer)[i], "NA", new.nt[[1]])
        # add new max.depth nt directly after the originial nt
        original_id <- which(names(urGrammar) == nt_double_pointer_df[i, 1])
        urGrammar <- c(urGrammar[1:original_id], new.nt, urGrammar[(original_id + 1):length(urGrammar)])
      }
    }
  }
  # initialise search table lenghts
  sear_dim <- data.frame("non-terminal" = names(urGrammar),
                         "rep" = 1,
                         stringsAsFactors = FALSE)
  #find all <> values and produce table with frequency
  sear_freq <- table(unlist(
    lapply(urGrammar,
           function(x) {
             # take unique values, only in case that a nt option has a nt occuring more than once take that number
             if(is.data.frame(x)) x <- as.character(x[1,])
             uGs <- unique(unlist(regmatches(x, gregexpr("<(.*?)>", x))))
             uGs_table <- table(uGs)
             uGa <- regmatches(x, gregexpr("<(.*?)>", x))
             uGa_table <- unlist(lapply(uGa, table))
             uGa_names <- names(uGa_table)
             for(name in unique(uGa_names)){
               uGs_table[names(uGs_table) == name] <- max(uGa_table[uGa_names == name],
                                                          uGs_table[names(uGs_table) == name])
             }
             return(rep(uGs, uGs_table))
           }
    )
  ))
  # update search table length with frequency
  tryCatch({
    for(i in seq_along(sear_freq)){
      ntfrq_id <- sear_dim$non.terminal == names(sear_freq)[i]
      sear_dim$rep[ntfrq_id] <- sear_freq[i]
    }
  }, error = function(err) {
    stop("Error: Number of defined non-terminals and number of used non-terminals in grammar is not the same! Check if there are any typos in your Grammar",
         call. = FALSE)
  })
  # check if there is one refered nt in the grammar that has less repetition then the one refering -> update
  for(i in seq_along(urGrammar)){
    current_nt <- names(urGrammar)[i]
    refered_nts <- unique(unlist(regmatches(as.character(urGrammar[[i]]), gregexpr("<(.*?)>", urGrammar[[i]]))))
    for(k in seq_along(refered_nts)){
      value_current <- sear_dim$rep[sear_dim$non.terminal == current_nt]
      value_refered <- sear_dim$rep[sear_dim$non.terminal == refered_nts[k]]
      if(value_current > value_refered) {
        sear_dim$rep[sear_dim$non.terminal == refered_nts[k]] <- value_current
      }
    }
    
  }
  # recursive repetitions of non-terminals need to have the same amount of reps as the original
  if(length(r_nt) > 0){
    # for double recursion
    if(length(r_nt_double_rec) > 0){
      r_nt_double_rec <- do.call(rbind, r_nt_double_rec)
      r_nt_double_rec <- data.frame(nt = r_nt_double_rec[, 3], rep = as.numeric(r_nt_double_rec[, 2]),
                                    stringsAsFactors = FALSE)
      r_nt_double_rec <- r_nt_double_rec[r_nt_double_rec$nt == unique(r_nt_double_rec$nt), ]
      changed_sear_dim <- sear_dim$non.terminal %in% r_nt_double_rec$nt
      sear_dim$rep[changed_sear_dim] <-  sear_dim$rep[changed_sear_dim] - r_nt_double_rec$rep
    }
    
    for (i in seq_along(r_nt)){
      nt_chain_ind <- grep(substr(r_nt[i], 1, nchar(r_nt[i]) - 1) , sear_dim$non.terminal)
      sear_dim$rep[nt_chain_ind] <- max(sear_dim$rep[nt_chain_ind])
      
    }
    # add double recursion values again
    if(length(r_nt_double_rec) > 0){
      sear_dim$rep[changed_sear_dim] <-  sear_dim$rep[changed_sear_dim] + r_nt_double_rec$rep
    }
  }
  # double pointing repetitions of non-terminals need to have the same amount of reps as the original
  if(length(nt_double_pointer) > 0){
    for(i in 1:length(nt_double_pointer)){
      end_point <- names(nt_double_pointer)[i]
      maxdepth_start <- nt_double_pointer[[i]][, 3]
      new_rep <- sear_dim$rep[sear_dim$non.terminal == nt_double_pointer[[i]][, 1]]
      # in case that a nt is occuring more often in a single option, new_rep can be larger
      new_rep <- max(new_rep, sear_dim$rep[sear_dim$non.terminal %in% c(end_point, maxdepth_start)])
      sear_dim$rep[sear_dim$non.terminal %in% c(end_point, maxdepth_start)] <- new_rep
    }
  }
  # unfold table
  search_table <- data.frame("non-terminal" = rep(sear_dim$non.terminal, sear_dim$rep),
                             "value" = NA,
                             stringsAsFactors = FALSE)
  # in case that max.depth nts have only one option -> remove them, otherwise let them be
  
  # get num index of nts
  nr_nt <- gsub("[^0-9]", "", search_table$non.terminal)
  # which are max.depth nts
  rm_nt <- which(as.numeric(nr_nt) > max.depth)
  # Which terminals have only one option
  option_l1 <- logical(length(rm_nt))
  for(i in seq_along(rm_nt)){
    option_l1[i] <- length(unlist(urGrammar[search_table$non.terminal[rm_nt[i]]])) == 1
  }
  # remove those with only one option from the search table
  rm_nt <- rm_nt[option_l1]
  if(length(rm_nt) > 0){
    search_table <- search_table[-rm_nt, ]
  }
  row.names(search_table) <- NULL
  return(search_table)
}

apply.search <- function(Grammar, search_table, max.depth = NULL, modulo = FALSE, na_remove = TRUE){
  ### create function from a search_table vector ###
  
  # subset s4 object
  if(!is.list(Grammar)){
    Grammar <- Grammar@grammar
  }
  # remove space from grammar
  Grammar <- lapply(Grammar, function(x) gsub(" ", "", x))
  # find out max.depth from given search table
  tryCatch({
    nr_nt <- gsub("[^0-9]","", search_table$non.terminal)
  }, error = function(cond) {
    stop("The search table does not have an non.terminal column. Use the search.table function to generate a template table for your Grammar.",
         call. = FALSE)
  })
  # set max.depth only if there are recursive elements
  if(is.null(max.depth)){
    if(sum(nr_nt != "") > 0){
      max.depth <- max(as.numeric(nr_nt), na.rm = TRUE)
    } else max.depth <- 0
  }
  # Check if given search_table has the same diminsion as given by the search.table function
  synth_table <- search.table(Grammar, max.depth = max.depth)
  if(sum(dim(synth_table) != dim(search_table)) > 0){
    stop(paste0("The applied search table has the wrong dimensions. Try to specify the correct max.depth from the search.table function.",
                "In case this error still occures, check if the given search table was really produced by search.table()."))
  }
  # unroll Grammar
  urGrammar <- grammar.unroller(Grammar, max.depth = max.depth)
  # Recursive and double recusive nts 
  # Recursive or not?
  nonterm <- names(Grammar)
  rec <- vector(mode ="list", length = length(Grammar))
  for(i in 1:length(Grammar)){
    rec[[i]] <- grep(nonterm[i], Grammar[[i]])
  }
  # Which are the recursive non-terminals
  r_nt <- names(Grammar)[sapply(rec, function(x) length(x) > 0)]
  # Do the recursive non-terminals include other recursive non-terminals?
  # loop over all recursive non-terminals to investigate if they have other rec. nts
  r_nt_pointer <- vector(mode = "list", length = length(r_nt))
  names(r_nt_pointer) <- r_nt
  r_nt_double_rec <- vector(mode = "list", length = length(r_nt))
  names(r_nt_double_rec) <- r_nt
  # define Grammar hierachy
  gHierachy <- data.frame(nt = names(Grammar), hierachy = 1:length(Grammar))
  number_of_recursives <- length(r_nt_pointer)
  if(number_of_recursives > 0){
    for(l in seq_along(r_nt)){
      # loop over all recursive non-terminals that are not investigated
      for(n in seq_along(r_nt[-l])) {
        is_rec.nt_inrec.nt <- grepl(r_nt[-l][n], Grammar[r_nt[l]])
        if(is_rec.nt_inrec.nt) {r_nt_pointer[[l]] <- c(r_nt_pointer[[l]], r_nt[-l][n])}
      }
    }
    # remove non-terminals that have no other recursive non-terminals included in their statement
    r_nt_pointer <- r_nt_pointer[!unlist(lapply(r_nt_pointer, is.null))]
    
  }
  number_of_recursives <- length(r_nt_pointer)
  if(number_of_recursives > 0){
    for (k in 1:number_of_recursives){
      NTi <- r_nt_pointer[[k]]
      NTi_name <- names(r_nt_pointer)[k]
      # is one of the pointed recursive elements lower in the hierachy (later in Grammar)?
      hierachy_of_nt <- gHierachy[gHierachy$nt == NTi_name, "hierachy"]
      hierachy_of_pointer <- gHierachy[gHierachy$nt %in% NTi, "hierachy"]
      r_nt_double_rec[[k]] <- NTi[hierachy_of_nt > hierachy_of_pointer]
    }
    # remove non-terminals that have no double recursion
    r_nt_double_rec <- r_nt_double_rec[unlist(lapply(r_nt_double_rec, function(x) length(x) > 0))]
    # Do only if there is a double recursion present
    if(length(r_nt_double_rec) > 0){
      # change double recursion links: all that point upwards to a recursion are set to max.depth
      for(i in 1:length(r_nt_double_rec)){
        nt_to_change <- names(r_nt_double_rec)[i]
        nt_to_change_into <- r_nt_double_rec[[i]]
        nts_to_change <- names(urGrammar)[grep(nt_to_change, gsub("[1-9]", "", names(urGrammar)))]
        for(k in seq_along(nts_to_change)){
          new_nt_name <- paste0(substr(nt_to_change_into, 1, nchar(nt_to_change_into) - 1), max.depth + 1, ">")
          urGrammar[[nts_to_change[k]]] <- gsub(nt_to_change_into, new_nt_name, urGrammar[[nts_to_change[k]]])
        }
      }
    }
  }
  # Double pointing without recursion 
  #   Find out if there are any nts that point at each other and create a infinite recursion.
  #   Don't take in those that are already dealt with in the double recursive part.
  nt_pointer <- vector(mode = "list", length = length(Grammar))
  names(nt_pointer) <- nonterm
  if(length(r_nt_pointer) > 0) {
    nt_pointer <- nt_pointer[!(names(nt_pointer) %in% r_nt_pointer)]
  }
  pointer <- names(nt_pointer)
  for(l in seq_along(nt_pointer)){
    # loop over all recursive non-terminals that are not investigated
    for(n in seq_along(pointer[-l])) {
      is_rec.nt_inrec.nt <- grepl(pointer[-l][n], Grammar[pointer[l]])
      if(is_rec.nt_inrec.nt) nt_pointer[[l]] <- c(nt_pointer[[l]], pointer[-l][n])
    }
  }
  # in case there is a double pointer, the higher non-terminal will be set to the value at max.depth
  # Do only if double pointer exist
  nt_pointer <- nt_pointer[!unlist(lapply(nt_pointer, is.null))]
  nt_double_pointer <- vector(mode = "list", length(nt_pointer))
  names(nt_double_pointer) <- names(nt_pointer)
  for (i in seq_along(nt_pointer)){
    NTi <- nt_pointer[[i]]
    NTi_name <- names(nt_pointer)[i]
    # is one of the pointed recursive elements lower in the hierachy (later in Grammar)?
    hierachy_of_nt <- gHierachy[gHierachy$nt == NTi_name, "hierachy"]
    hierachy_of_pointer <- gHierachy[gHierachy$nt %in% NTi, "hierachy"]
    nt_double_pointer[[i]] <- NTi[hierachy_of_nt > hierachy_of_pointer]
  }
  # remove non-terminals that are higher in hierachy that the ones they point to
  nt_double_pointer <- nt_double_pointer[unlist(lapply(nt_double_pointer, function(x) length(x) > 0))]
  # change double recursion links: all that point upwards to a recursion are set to max.depth
  if(length(nt_double_pointer) > 0){
    number_of_changes <- numeric()
    for(i in 1:length(nt_double_pointer)){
      nt_to_change <- names(nt_double_pointer)[i]
      nt_to_change_into <- nt_double_pointer[[i]]
      nts_to_change <- names(urGrammar)[grep(nt_to_change, gsub("[1-9]", "", names(urGrammar)))]
      # get highest number of current Grammar for nts_to_change_into
      search_nt_to_change_into <- substr(nt_to_change_into, 1, nchar(nt_to_change_into) - 1)
      max.change <- names(urGrammar)[grep(search_nt_to_change_into, names(urGrammar))]
      max.change <- gsub(search_nt_to_change_into, "", max.change)
      max.change <- gsub(">", "", max.change)
      if(sum(max.change %in% "") != length(max.change)){
        depth.id <- max(as.numeric(max.change), na.rm = TRUE)
      } else depth.id <- 1
      
      for(k in seq_along(nts_to_change)){
        new_nt_name <- paste0(substr(nt_to_change_into, 1, nchar(nt_to_change_into) - 1), depth.id, ">")
        urGrammar[[nts_to_change[k]]] <- gsub(nt_to_change_into, new_nt_name, urGrammar[[nts_to_change[k]]])
      }
      number_of_changes <- c(number_of_changes, length(nts_to_change))
    }
    # produce list with added number of max.depth objects due to double recursion
    for(i in 1:length(nt_double_pointer)) {
      nt_double_pointer[[i]] <- cbind(nt_double_pointer[[i]], number_of_changes[i])
      nt_double_pointer[[i]] <- cbind(nt_double_pointer[[i]],
                                      paste0(substr(nt_double_pointer[[i]][1], 1,
                                                    nchar(nt_double_pointer[[i]][1]) - 1),
                                             depth.id,
                                             ">"))
    }
    # Add new terms max.depth terms to Grammar
    if(length(nt_double_pointer) > 0){
      nt_double_pointer_df <- do.call(rbind, nt_double_pointer)
      for(i in 1:nrow(nt_double_pointer_df)){
        new.nt <- vector(mode = "list", length = 1)
        names(new.nt) <- nt_double_pointer_df[i, 3]
        new.nt[[1]] <- urGrammar[[nt_double_pointer_df[i, 1]]]
        # remove all nts
        new.nt[[1]] <- gsub("<.*?>", "NA", new.nt[[1]])
        # add new max.depth nt directly after the originial nt
        original_id <- which(names(urGrammar) == nt_double_pointer_df[i, 1])
        urGrammar <- c(urGrammar[1:original_id], new.nt, urGrammar[(original_id + 1):length(urGrammar)])
      }
    }
  }
  
  #create list from search table
  tab <- split(search_table, search_table$non.terminal)
  # reduce list to only include integer vectors
  tab <- lapply(tab, function(x) x[,2])
  # compare number of possible options and the max chosen option for each non-terminal
  # in case the max is larger than the possible number of options -> error
  args_in_tab <- unlist(lapply(tab, max))
  if(sum(is.na(args_in_tab)) > 0) stop("Please provide a search table without NA values.")
  args_in_grammar <- unlist(lapply(urGrammar, length))
  args_in_grammar <- args_in_grammar[names(args_in_grammar) %in% names(args_in_tab)]
  args_in_tab <- data.frame(name = names(args_in_tab), rep = args_in_tab)
  args_in_grammar <- data.frame(name = names(args_in_grammar), rep = args_in_grammar)
  args_grammar_tab <- merge(args_in_grammar, args_in_tab, by = "name", suffixes = c(".grammar", ".tab"))
  if(modulo == FALSE){
    if(sum(args_grammar_tab$rep.grammar < args_grammar_tab$rep.tab) > 0){
      stop("There is at least one value in the search table, that is
         larger then the numbers of options for that non-terminal!\n
         e.g. Grammar:      <a> = <b> | <c>
              search table: <a> = 3")
    }
  } else {
    # goe through every tab value, if value is larger the max options use the modulo
    args_in_grammar <- args_grammar_tab$rep.grammar
    names(args_in_grammar) <- args_grammar_tab$name
    for(i in 1:length(tab)){
      for(j in seq_along(tab[[i]])){
        if(tab[[i]][j] > args_in_grammar[names(tab)[i]]){
          tab[[i]][j] <- tab[[i]][j] %% args_in_grammar[names(tab)[i]]
          if(tab[[i]][j] == 0) tab[[i]][j] <- args_in_grammar[names(tab)[i]]
        }
      }
    }
  }
  # Initial Grammar subsetting
  rgram <- as.character(urGrammar[[1]][tab[[nonterm[1]]][1]])
  # delete first tab value
  tab[[nonterm[1]]]<- tab[[nonterm[1]]][-1]
  # While loop until rgram does not includes anymore non-terminal variables
  while(grepl("<(.*?)>", rgram)){
    # get non-terminals in current rgram
    nt_in_rgram <- unlist(regmatches(rgram, gregexpr("<(.*?)>", rgram)))
    # find out if any are at the max.depth length
    terminals_ind <- which(as.numeric(gsub("[^0-9]","",nt_in_rgram)) > max.depth)
    # if there are terminals, split them from the non-terminals
    if(length(terminals_ind) > 0){
      t_in_rgram <- nt_in_rgram[terminals_ind]
      nt_in_rgram <- nt_in_rgram[-terminals_ind]
      for(j in seq_along(t_in_rgram)){
        # if the terminal is in urGrammar than take one of those and put NAs instead of nts
        if(t_in_rgram[j] %in% names(tab)){
          # get Grammar value chosen in tab
          t_value <- urGrammar[[t_in_rgram[j]]][tab[[t_in_rgram[j]]][1]]
          # delete used tab value
          tab[[t_in_rgram[j]]] <- tab[[t_in_rgram[j]]][-1]
          # change all remaining nts to NA
          t_value <- gsub("<.*?>", "NA", t_value)
          # put it in
          rgram <- gsub(t_in_rgram[j], t_value, rgram)
        } else {
          # get Grammar value
          t_value <- urGrammar[[t_in_rgram[j]]]
          # change all remaining nts to NA
          t_value <- gsub("<.*?>", "NA", t_value)
          # put it in
          rgram <- gsub(t_in_rgram[j], t_value, rgram)
        }
      }
      if(length(nt_in_rgram) == 0) next
    }
    # get replacements for the non-terminals
    # tab value:
    tab_val <- lapply(tab[nt_in_rgram], function(x) x[1])
    # delete already used tab values -> onyl the first will ever be used
    tab[nt_in_rgram] <- lapply(tab[nt_in_rgram], function(x) x[-1])
    # add non-terminal name to value list -> makes next step easier
    for(i in 1:length(tab_val)) tab_val[[i]] <- data.frame(nt = names(tab_val)[i], val = tab_val[[i]],
                                                           stringsAsFactors = FALSE)
    # get replacement
    nt_replacement <- lapply(tab_val, function(x) urGrammar[[x$nt]][x$val])
    # give replacement terminals a unique name, so that they are changed only once
    for(j in 1:length(nt_replacement)) nt_replacement[[j]] <- gsub(">", "*>", nt_replacement[[j]])
    # replace non-terminals with replacement
    for(i in seq_along(nt_in_rgram)){
      rgram <- sub(nt_in_rgram[i], nt_replacement[[i]], rgram)
    }
    # remove unique markers
    rgram <- gsub("\\*>", ">", rgram)
    # end of while loop
  }
  if (na_remove) rgram <- na.remover(rgram)
  return(rgram)
}

simplify_functions <- function(fun, function_variables){
  ### simplifies a mathematical expression given as a string ###
  
  if (length(fun) > 1) stop("fun should be one character string. If you want to apply the function to a vector of fun use parlapply_simplify")
  # remove spaces
  fun <- gsub(" ", "", fun)
  nums <- length(grep("numeric", fun)) > 0 # are there "numeric" values
  if(nums){
    # If there are no function_variables -> return NA
    tf_splitted <- unlist(strsplit(fun, "(?=[+-/*)()])", perl = TRUE))
    if(sum(function_variables %in% tf_splitted) == 0) return(NA)
    # here starts the real function
    # which numerics are only sourrounded by +/- -> reduce it to just one + numeric
    tf_len <- nchar(fun)
    while(TRUE){
      now_len <- nchar(fun)
      if( length(grep("+numeric+", fun, fixed = TRUE)) > 0) fun <- gsub("+numeric+", "+", fun, fixed = TRUE)
      if( length(grep("-numeric+", fun, fixed = TRUE)) > 0) fun <- gsub("-numeric+", "+", fun, fixed = TRUE)
      if( length(grep("+numeric-", fun, fixed = TRUE)) > 0) fun <- gsub("+numeric-", "-", fun, fixed = TRUE)
      if( length(grep("-numeric-", fun, fixed = TRUE)) > 0) fun <- gsub("-numeric+", "+", fun, fixed = TRUE)
      
      if(substr(fun, nchar(fun) - 7, nchar(fun)) == "+numeric") {
        fun <- fun <- substr(fun, 1, nchar(fun) - 8)
      }
      if(substr(fun, nchar(fun) - 7, nchar(fun)) == "-numeric") {
        fun <- fun <- substr(fun, 1, nchar(fun) - 8)
      }
      if(substr(fun, 1, 8) == "numeric+") {
        fun <- substring(fun, 9)
      }
      if(substr(fun, 1, 8) == "numeric-") {
        fun <- substring(fun, 9)
      }
      if(now_len == nchar(fun)) break
    }
    # in case some numerics were removed -> add one +numeric
    if(tf_len != nchar(fun)) fun <- paste0(fun, "+numeric")
    
    for(i in 1:4){
      # replace obvious numerics
      fun <- gsub("numeric/numeric", "numeric", fun, fixed = TRUE)
      fun <- gsub("numeric*numeric", "numeric", fun, fixed = TRUE)
      
      # if a single numeric is between parentheses -> remove parantheses
      fun <- gsub("(numeric)", "numeric", fun, fixed = TRUE)
    }
    # split fun at +-
    fun <- unlist(strsplit(fun, "(?=[+-])", perl = TRUE))
    fun <- fun[fun != "0"]
    # are two numerics seperated by a /
    for(i in 1:length(fun)){
      sfun <- unlist(strsplit(fun[i], "(?=[/])", perl = TRUE))
      if(length(sfun) > 1){
        for(k in 1:length(sfun)) sfun[k] <- gsub("numeric", paste0("numeric", k), sfun[k])
        fun[i] <- paste0(sfun, collapse = "")
      }
    }
    for(i in 1:length(fun)) fun[i] <- gsub("numeric", paste0("numeric", i), fun[i])
    fun <- paste0(fun, collapse = "")
    if(substr(fun, 1, 1) %in% c("+", "-")) fun <- substring(fun, 2)
    if(substring(fun, nchar(fun)) %in% c("+", "-")) fun <- substr(fun, 1, nchar(fun)-1)
  }
  fun <- Deriv::Simplify(fun)
  fun <- Deriv::Simplify(fun)
  fun <- gsub(" ", "", fun)
  if(nums){
    fun <- gsub("numeric+[0-9]", "numeric", fun)
    fun <- gsub("numeric+[0-9]", "numeric", fun)
    fun <- gsub("numeric+[0-9]", "numeric", fun)
    # change numeric with power values to numerics
    fun <- gsub("numeric+\\^+[0-9]", "numeric", fun)
    fun <- gsub("numeric/numeric", "numeric", fun, fixed = TRUE)
    fun <- gsub("numeric*numeric", "numeric", fun, fixed = TRUE)
  }
  # clean numerics one more time
  if(nums){
    # Check if there are any obviously non important nums
    while(TRUE){
      now_len <- nchar(fun)
      if( length(grep("(numeric+numeric)", fun, fixed = TRUE)) > 0) fun <- gsub("numeric+numeric", "numeric", fun, fixed = TRUE)
      if( length(grep("(numeric-numeric)", fun, fixed = TRUE)) > 0) fun <- gsub("numeric+numeric", "numeric", fun, fixed = TRUE)
      if(now_len == nchar(fun)) break
    }
    tf_len <- nchar(fun)
    while(TRUE){
      now_len <- nchar(fun)
      if( length(grep("+numeric+", fun, fixed = TRUE)) > 0) fun <- gsub("+numeric+", "+", fun, fixed = TRUE)
      if( length(grep("-numeric+", fun, fixed = TRUE)) > 0) fun <- gsub("-numeric+", "+", fun, fixed = TRUE)
      if( length(grep("+numeric-", fun, fixed = TRUE)) > 0) fun <- gsub("+numeric-", "-", fun, fixed = TRUE)
      if( length(grep("-numeric-", fun, fixed = TRUE)) > 0) fun <- gsub("-numeric+", "+", fun, fixed = TRUE)
      if(substr(fun, nchar(fun) - 7, nchar(fun)) == "+numeric") {
        ffun <- substr(fun, 1, nchar(fun) - 8)
      }
      if(substr(fun, nchar(fun) - 7, nchar(fun)) == "-numeric") {
        fun <- substr(fun, 1, nchar(fun) - 8)
      }
      if(substr(fun, 1, 8) == "numeric+") {
        fun <- substring(fun, 9)
      }
      if(substr(fun, 1, 8) == "numeric-") {
        fun <- substring  (fun, 9)
      }
      if(now_len == nchar(fun)) break
    }
    # in case some numerics were removed -> add one +numeric
    if(tf_len != nchar(fun)) fun <- paste0(fun, "+numeric")
  }
  # did simplifying remove spatial predictors? test if there are still some
  # If there are no spatial_predictors -> return NA
  tf_splitted <- unlist(strsplit(fun, "(?=[+-/*()^])", perl = TRUE))
  if(sum(function_variables %in% tf_splitted) == 0) return(NA)
  fun <- gsub("*", " * ", fun, fixed = TRUE)
  fun <- gsub("+", " + ", fun, fixed = TRUE)
  fun <- gsub("-", " - ", fun, fixed = TRUE)
  return(fun)
}

parlapply_simplify <- function(funs, function_variables, no_cores = NULL){
  ## parallel applying simplify_functions on a vector of functions ###
  
  # Calculate the number of cores
  if(is.null(no_cores)){
    no_cores <- parallel::detectCores() - 2
  }
  # Initiate cluster
  start <- Sys.time()
  functions_list <- parallel::mclapply(funs,
                                       simplify_functions,
                                       function_variables = function_variables,
                                       mc.cores = no_cores)
  end <- Sys.time()
  cat("\nParallel simplify took ", (as.numeric(end) - as.numeric(start))/60, "minutes\n")
  return(unlist(functions_list))
}

rand.grammar.sampler <- function(gram, gram_table, depth){
  ## samples a random function from a grammar ###
  
  lens <- unlist(lapply(gram, function(x) length(x)))
  names(lens) <- gsub("[<, >]", "", names(lens))
  for(i in 1:nrow(gram_table)){
    nt_index <- gsub("[<, >, 0-9]", "", gram_table[i, 1]) == names(lens)
    gram_table[i, 2] <- sample(lens[nt_index], 1)
  }
  return(list(gram_vector = gram_table$value,
              "transfer_function" = apply.search(gram, gram_table, max.depth = depth)
  ))
}

functions.creator <- function(n, gram, gram_table, depth){
  ### samples n random functions from a grammar ###
  
  functions <- vector(mode = "list", length = n)
  for(i in 1:n){
    functions[[i]] <- rand.grammar.sampler(gram = gram, gram_table = gram_table, depth = depth)
  }
  gram_vectors <- matrix(NA, nrow = nrow(gram_table), ncol = n)
  functions_vector <- character(n)
  for(i in 1:n){
    gram_vectors[, i] <- functions[[i]]$gram_vector
    functions_vector[i] <- functions[[i]]$transfer_function
  }
  
  return(list("gram_vectors" = gram_vectors,
              "functions_vector" = functions_vector))
}

par.grammar.sampler <- function(cfgram, max.depth, n, save_feather = TRUE, parallel = TRUE, no_cores = NULL){
  ### parallel sampling of random functions from a grammar ###
  
  # prepare grammar and table
  gram_table <- search.table(cfgram, max.depth = max.depth)
  gram <- cfgram@grammar
  if(parallel){
    # Calculate the number of cores
    if(is.null(no_cores)){
      no_cores <- parallel::detectCores() - 2
    }
    ns <- rep(as.integer(n/no_cores), no_cores)
    # Parallel run with mclapply
    start <- Sys.time()
    cat("Run startet at: ", as.character(start))
    RNGkind("L'Ecuyer-CMRG")
    functions_list <- parallel::mclapply(ns, functions.creator, mc.cores = no_cores,
                               gram = gram,
                               gram_table = gram_table,
                               depth = max.depth,
                               mc.set.seed = TRUE)
    end <- Sys.time()
    cat("\nParallelized run with ", sum(ns), " iterations needed ", (end-start)/60/60, "hours.")
    #stopCluster(cl)
    time_name <- substr(as.character(end), 1, 10)
    gram_vectors <- functions_list[[1]]$gram_vector
    functions_vector <- functions_list[[1]]$functions_vector
    for(i in 2:no_cores){
      gram_vectors <- cbind(gram_vectors, functions_list[[i]]$gram_vector)
      functions_vector <- c(functions_vector, functions_list[[i]]$functions_vector)
    }
  } else {
    start <- Sys.time()
    cat("Run startet at: ", as.character(start))
    functions <- functions.creator(n, gram = gram, gram_table = gram_table, depth = max.depth)
    end <- Sys.time()
    cat("\nNot parallelized run with", n, "iterations needed", (as.numeric(end)-as.numeric(start))/60/60, "hours.")
    time_name <- substr(as.character(end), 1, 10)
    gram_vectors <- functions$gram_vectors
    functions_vector <- functions$functions_vector
    
  }
  functions_list <- list("Grammatic_vectors" = gram_vectors,
                         "transfer_functions" = functions_vector)
  # save as feather
  if(save_feather){
    feather::write_feather(functions_list, paste0("function_space", "-", time_name, ".csv"), row.names = FALSE)
  }
  functions_df <- data.frame(t(functions_list$Grammatic_vectors), functions_list$transfer_functions,
                             stringsAsFactors = FALSE)
  names(functions_df) <- c(gram_table$non.terminal, "Transfer_Function")
  return(functions_df)
}
