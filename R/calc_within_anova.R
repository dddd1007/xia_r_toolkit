library(magrittr)

calc_within_anova <- function(filtered_data, subject_num, factors, dep = "RT") {
    # 确保 factors 包括 subject_num
    grouping_columns <- c(subject_num, factors)

    # 排除三个标准差以外的数据
    filtered_data <- filtered_data %>%
        dplyr::group_by_at(grouping_columns) %>%
        dplyr::filter(abs(get(dep) - mean(get(dep), na.rm = TRUE)) < (sd(get(dep), na.rm = TRUE) * 3))

    # 按照指定的分组计算均值
    mean_data <- filtered_data %>%
        dplyr::group_by_at(grouping_columns) %>%
        dplyr::summarise(mean_rt = mean(get(dep), na.rm = TRUE)) %>%
        dplyr::as_tibble()

    # 进行方差分析
    formula <- as.formula(paste("mean_rt ~ factor(", paste(factors, collapse = ") * factor("), ")", sep = ""))
    mod <- rstatix::Anova(lm(formula, data = mean_data), type = 3)

    # 打印 ANOVA 结果并返回
    print(mod)
    return(mod)
}
