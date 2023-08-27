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
        anova_test(
            dv = mean_rt, wid = within, # nolint
            within = factors,
            type = type
        ) -> rstatix_anova
    get_anova_table(rstatix_anova) -> rstatix_results

    # 进行事后检验
    post_hoc <- tukey_hsd(mean_data, formula)

    ### 方法3 使用 ezANOVA
    if (debug) {
        print(head(input_data))
    }

    if (sum(colnames(input_data) == dep) == 0 &
        sum(factors %in% colnames(input_data)) != length(factors)) {
        stop("dep or factors is not in input_data")
    } else {
        ez_results <- rlang::eval_tidy(rlang::expr(ez::ezANOVA(
            data = input_data,
            dv = .(!!sym(dep)),
            wid = .(!!sym(within)),
            within = .(!!!rlang::syms(factors)),
            detailed = TRUE
        )))
    }

    if (print_desc_report) {
        cat("\033[31m=== Descibtion ===\033[39m\n")
        input_data %>%
            group_by(dplyr::across(all_of(factors))) %>%
            report::report() %>%
            summary() %>%
            print()
        cat("\033[31m=== End Report ===\033[39m\n")
    }

    if (generate_report) {
        apaTables::apa.aov.table(lm_output, filename = paste0(report_path, "/anova_table.doc"))
        apaTables::apa.ezANOVA.table(ez_results, filename = paste0(report_path, "/ezANOVA_table.doc"))
        print(paste("The anova report are saved to", getwd()))
    }

    return(list(
        "car_results" = car_results,
        "rstatix_results" = rstatix_results,
        "ez_results" = ez_results,
        "post_hoc" = post_hoc
    ))
}

#' pc_simon_plot
#'
#' 这个函数用于绘制 Simon 任务的反应时间数据的条形图。
#'
#' @param raw_data 原始数据，应为一个数据框
#' @param dep_var 依赖变量，应为一个字符串，表示数据框中的一列
#' @param factors 因子，应为一个字符向量，表示数据框中的多列
#' @return 返回一个 ggplot 对象，表示绘制的条形图
#' @export
#' @examples
#' # 假设我们有一个名为 df 的数据框，其中包含因子 "subject_num"、"condition" 和 "trial"，以及依赖变量 "RT"
#' pc_simon_plot(df, dep_var = "RT", factors = c("subject_num", "condition", "trial"))
pc_simon_plot <- function(input_data, dep_var, factors,
                          facet_by = "condition",
                          color_by = "congruency",
                          ylim,
                          debug = FALSE) {
    library(tidyr)

    anova_table <- input_data %>%
        group_by(!!!syms(factors)) %>%
        dplyr::summarise(mean_dep_var = mean(!!sym(dep_var))) %>%
        as_tibble()

    data_forplot <- unite(anova_table,
        col = "combined_factor",
        !!!syms(setdiff(factors, facet_by)), sep = "/", remove = FALSE
    )

    data_forplot$combined_factor <- as.factor(data_forplot$combined_factor)

    if (debug) {
        # 将 tibble data_forplot 转换为 dataframe
        data_forplot <- as.data.frame(data_forplot)
        print(data_forplot)
    }

    anova_plot_bar <- ggplot(data_forplot, aes(x = combined_factor, y = mean_dep_var, fill = !!sym(color_by))) + # nolint
        geom_bar(stat = "identity", position = "dodge") +
        facet_grid(paste0(". ~ ", facet_by)) +
        geom_errorbar(
            aes(
                ymin = mean_dep_var - sd(mean_dep_var),
                ymax = mean_dep_var + sd(mean_dep_var)
            ),
            width = 0.2,
            color = "black"
        ) +
        labs(y = paste0(dep_var, " (ms)")) +
        ggpubr::theme_pubr() +
        theme(text = element_text(size = 14)) +
        coord_cartesian(ylim = ylim)

    return(anova_plot_bar)
}
