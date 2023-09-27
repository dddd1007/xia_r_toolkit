library(tidyverse)
library(magrittr)
library(rstatix)
library(ggplot2)

#' calc_within_anova
#'
#' 计算组内方差分析
#'
#' @param raw_data 原始数据，应为一个数据框
#' @param within_factor 组内因子，默认为 "subject_num"
#' @param factors 其他因子，应为一个字符向量
#' @param dep 依赖变量，默认为 "RT"
#' @param type 决定 Anova 的计算类别
#' @param remove_outlier 是否移除离群值，默认为 TRUE
#' @return 返回一个包含 ANOVA 结果的对象
#' @export
#' @examples
#' # 假设我们有一个名为 df 的数据框，其中包含因子 "subject_num" 和 "condition"，以及依赖变量 "RT"
#' calc_within_anova(df, within_factor = "subject_num", factors = c("condition"), dep = "RT")
calc_within_anova <- function(input_data,
                              within = "subject_num",
                              factors,
                              dep = "RT",
                              type = 3,
                              generate_report = FALSE,
                              report_path = NULL,
                              print_desc_report = TRUE,
                              debug = FALSE) {
    grouping_columns <- c(within, factors)
    if (debug) {
        print(head(input_data))
    }

    ### 方法1 使用car包进行分析
    # 按照指定的分组计算均值
    mean_data <- input_data %>%
        dplyr::group_by_at(grouping_columns) %>%
        dplyr::summarise(mean_rt = mean(get(dep), na.rm = TRUE)) %>%
        dplyr::as_tibble()

    # 将 mean_data 中列名在 factors变量中的列都转换为 factor 类型
    mean_data[factors] <- lapply(mean_data[factors], factor)

    # 进行方差分析
    formula <- as.formula(paste("mean_rt ~ ", paste(factors, collapse = " * "), sep = ""))
    lm_output <- lm(formula, data = mean_data)
    car_results <- car::Anova(lm_output, type = type)

    ### 方法2 使用 rstarix 进行方差分析
    mean_data %>%
        rstatix::anova_test(
            dv = mean_rt, wid = within, # nolint
            within = factors,
            type = type,
            effect.size = "pes"
        ) -> rstatix_anova
    rstatix::get_anova_table(rstatix_anova) -> rstatix_results

    # 进行事后检验
    post_hoc <- tukey_hsd(mean_data, formula)

    ### 方法3 使用 ezANOVA
    if (sum(colnames(input_data) == dep) == 0 ||
        sum(factors %in% colnames(input_data)) != length(factors)) {
        stop("dep or factors is not in input_data")
    } else {
        # Use ez::ezANOVA to calculate ANOVA and get partial eta square
        ez_results <- rlang::eval_tidy(rlang::expr(ez::ezANOVA(
            data = input_data,
            dv = .(!!sym(dep)),
            wid = .(!!sym(within)),
            within = .(!!!rlang::syms(factors)),
            detailed = TRUE,
            type = type
        )))
        # Extract partial eta square from ez_results
        ez_results$ANOVA$"partial_eta_square" <- ez_results$ANOVA$SSn / (ez_results$ANOVA$SSn + ez_results$ANOVA$SSd)
    }

    if (print_desc_report) {
        cat("\033[31m=== 开始报告统计结果 ===\033[39m\n")
        cat("\033[34m=== 数据描述统计 :\033[39m\n")
        mean_data %>%
            group_by(dplyr::across(all_of(factors))) %>%
            summarise(mean_rt = mean(mean_rt)) %>%
            print()
        cat("\033[31m=== 统计报告结束 ===\033[39m\n")
    }

    if (generate_report) {
        apaTables::apa.aov.table(lm_output, filename = paste0(report_path, "/", dep, "_anova_table.doc"))
        write.csv(ez_results$ANOVA, paste0(report_path, "/", dep, "_ezANOVA_table.csv"))
        print(paste("The anova report are saved to", getwd()))
    }

    return(list(
        "car_results" = car_results,
        "rstatix_results" = rstatix_results,
        "ez_results" = ez_results,
        "post_hoc" = post_hoc
    ))
}
