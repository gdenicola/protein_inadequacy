# ==============================================================
# Donor-country helpers for the donor bootstrap (Script 6; R4.3, R4.6)
#
# The main analysis (Scripts 1 and 3) gives strata without their own nutriR
# distribution the shape of the single nearest donor country: similarity in
# country-mean GDD protein intake for protein, and in national energy supply
# before downstream waste for energy.
#
# The donor bootstrap in Script 6 instead lets each country take one of its
# DONOR_K nearest donor countries. donor_rank_shapes() returns, for every
# stratum that needs a donor, the shape it gets from its r-th nearest donor
# (r = 1, ..., DONOR_K), using the same within-donor steps as Step 4: same sex
# and age group; else nearest age group, same sex; else same age group,
# opposite sex; else nearest age group, opposite sex.
# ==============================================================

DONOR_K <- 5

donor_age_start <- function(a) suppressWarnings(as.numeric(sub("^([0-9.]+)-.*", "\\1", a)))

# Shape of one donor (rows of `shapes` for that donor) for target sex s, age ag
lookup_shape_in_donor <- function(sd, s, ag) {
  opp <- if (s == "Males") "Females" else "Males"
  nearest <- function(pool) pool[which.min(abs(donor_age_start(pool$age_group) - donor_age_start(ag))), , drop = FALSE]
  same <- sd[sd$sex == s, , drop = FALSE]
  oppd <- sd[sd$sex == opp, , drop = FALSE]
  if (any(same$age_group == ag)) {
    r <- same[same$age_group == ag, , drop = FALSE][1, ]
  } else if (nrow(same) > 0) {
    r <- nearest(same)
  } else if (any(oppd$age_group == ag)) {
    r <- oppd[oppd$age_group == ag, , drop = FALSE][1, ]
  } else if (nrow(oppd) > 0) {
    r <- nearest(oppd)
  } else {
    return(NULL)
  }
  tibble::tibble(sex = s, age_group = ag, match_iso3 = r$iso3, best_dist = r$best_dist, cv = r$cv)
}

# inp:   list(targets, shapes, metric, donors) saved by Script 1 or Script 3
#        targets: strata needing a donor (iso3, sex, age_group)
#        shapes:  nutriR shapes (iso3, sex, age_group, best_dist, cv)
#        metric:  named numeric vector (names = iso3) of the similarity metric
#        donors:  iso3 codes allowed as donors
# first: optional named vector (names = iso3) with the donor used in the main
#        analysis; it is kept as rank 1 and ranks 2..k follow by distance
# Returns one row per target stratum and donor rank: iso3, sex, age_group,
# donor_rank, match_iso3, best_dist, cv
donor_rank_shapes <- function(inp, k = DONOR_K, first = NULL) {
  shapes <- inp$shapes %>%
    dplyr::mutate(age_group = as.character(age_group)) %>%
    dplyr::filter(!is.na(cv), !is.na(best_dist)) %>%
    dplyr::distinct(iso3, sex, age_group, .keep_all = TRUE)
  metric <- inp$metric
  donors <- sort(intersect(intersect(inp$donors, unique(shapes$iso3)), names(metric)[!is.na(metric)]))
  targets <- inp$targets %>%
    dplyr::mutate(age_group = as.character(age_group)) %>%
    dplyr::distinct(iso3, sex, age_group)

  # shape of every donor for every target sex-age key
  keys <- dplyr::distinct(targets, sex, age_group)
  lut <- dplyr::bind_rows(lapply(donors, function(dn) {
    sd <- shapes[shapes$iso3 == dn, , drop = FALSE]
    dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(j) lookup_shape_in_donor(sd, keys$sex[j], keys$age_group[j])))
  }))

  # the k nearest donors of every target country
  nn <- dplyr::bind_rows(lapply(unique(targets$iso3), function(iso) {
    if (!iso %in% names(metric) || is.na(metric[[iso]])) return(NULL)
    cand <- setdiff(donors, iso)
    cand <- cand[order(abs(unname(metric[cand]) - metric[[iso]]), cand)]
    if (!is.null(first) && iso %in% names(first) && !is.na(first[[iso]])) {
      cand <- c(first[[iso]], setdiff(cand, first[[iso]]))
    }
    cand <- cand[seq_len(min(k, length(cand)))]
    tibble::tibble(iso3 = iso, donor_rank = seq_along(cand), match_iso3 = cand)
  }))

  merge(as.data.frame(targets), as.data.frame(nn), by = "iso3") %>%
    dplyr::inner_join(lut, by = c("match_iso3", "sex", "age_group")) %>%
    dplyr::select(iso3, sex, age_group, donor_rank, match_iso3, best_dist, cv) %>%
    tibble::as_tibble()
}
