library(magrittr)
#' clean_data
#'
#' 这个函数用于移除数据中的离群值。离群值被定义为在均值的三个标准差之外的值。
#'
#' @param raw_data 原始数据，应为一个数据框
#' @param target_col 目标列，应为一个字符串，表示数据框中的一列
#' @param filtered_by 用于过滤的列，应为一个字符向量，表示数据框中的多列
#' @param resp_cor_col 响应相关列，应为一个字符串，表示数据框中的一列
#' @param resp_cor_value 响应相关值，应为一个数值，表示响应相关列的值
#' @param debug 是否开启调试模式，默认为 FALSE。如果为 TRUE，将打印被移除的离群值的信息
#' @return 返回一个数据框，其中已经移除了离群值
#' @export
#' @examples
#' # 假设我们有一个名为 df 的数据框，其中包含因子 "subject_num" 和 "condition"，以及依赖变量 "RT"
#' clean_data(df, target_col = "RT", filtered_by = c("subject_num", "condition"), resp_cor_col = "condition", resp_cor_value = 1, debug = TRUE)
clean_data <- function(raw_data, target_col, filtered_by, resp_cor_col, resp_cor_value, debug = FALSE) {
    filtered_data <- raw_data %>%
        na.omit() %>%
        filter(!!rlang::sym(resp_cor_col) == resp_cor_value) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(filtered_by))) %>%
        dplyr::mutate(
            group_mean = mean(!!sym(target_col), na.rm = TRUE),
            group_sd = sd(!!sym(target_col), na.rm = TRUE)
        ) %>%
        dplyr::mutate(is_outlier = abs(!!sym(target_col) - group_mean) > (group_sd * 3)) %>%
        ungroup() %>%
        filter(!is_outlier) %>%
        select(-c(group_mean, group_sd, is_outlier))

    if (debug) {
        raw_data %>%
            na.omit() -> removed_na_data
        print(paste("The number of removed NA data is:", nrow(raw_data) - nrow(removed_na_data)))
        raw_data %>%
            na.omit() %>%
            filter(!!rlang::sym(resp_cor_col) == resp_cor_value) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(filtered_by))) %>%
            dplyr::mutate(
                group_mean = mean(!!sym(target_col), na.rm = TRUE),
                group_sd = sd(!!sym(target_col), na.rm = TRUE)
            ) %>%
            dplyr::mutate(is_outlier = abs(!!sym(target_col) - group_mean) > (group_sd * 3)) %>%
            dplyr::ungroup() %>%
            filter(is_outlier) -> droped_data
        print("=== The following data are dropped ===")
        print(head(droped_data))
        print(paste("The number of dropped data is:", nrow(droped_data)))
    }
    return(filtered_data)
}
