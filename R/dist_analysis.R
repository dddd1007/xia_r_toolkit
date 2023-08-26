#' delta_plot
#'
#' 计算 Simon 效应并绘制 delta 图
#'
#' @param input_data 输入数据，应为一个数据框
#' @param subject_col 被试列名，默认为 "subject_num"
#' @param congruency_col 一致性列名，默认为 "congruency"
#' @param rt_col 反应时间列名，默认为 "RT"
#' @param color_by 颜色分组列名
#' @param facet_by 分面分组列名
#' @param inc_value 一致性值，默认为 "inc"
#' @param factors_col 其他因子列名，应为一个字符向量
#' @param n_bins 分箱数量
#' @param remove_outlier 是否移除离群值，默认为 TRUE
#' @return 返回一个 ggplot 对象，表示 delta 图
#' @export
#' @examples
#' # 假设我们有一个名为 df 的数据框，其中包含因子 "subject_num" 和 "condition"，以及依赖变量 "RT"
#' delta_plot(df, subject_col = "subject_num", congruency_col = "congruency", rt_col = "RT", color_by = "condition", facet_by = "condition", factors_col = c("condition"), n_bins = 5)
delta_plot <- function(input_data, subject_col, congruency_col, rt_col,
                       color_by, facet_by,
                       inc_value = "inc", factors_col, n_bins, debug = FALSE) {
    library(tidyr)
    library(dplyr)

    # 获取 congruency_col 中不是 inc_value 的另一个值，如果有多于两个值则报错
    congruency_values <- unique(input_data[[congruency_col]])
    stopifnot(length(congruency_values) == 2, inc_value %in% congruency_values)
    con_value <- congruency_values[congruency_values != inc_value][1]

    # 对每个被试的数据先按照 factors 和 congruency_col 进行分组
    # 按照从小到大进行排序，然后平均分为 n_bins 个 bin
    # 每个 bin 的每个条件下计算平均值
    # Calculate global bin breaks

    input_data %>%
        dplyr::select(all_of(c(subject_col, factors_col, congruency_col, rt_col))) %>%
        dplyr::group_by(dplyr::across(all_of(c(factors_col, congruency_col, subject_col)))) %>%
        dplyr::arrange(!!sym(rt_col)) %>%
        dplyr::mutate(bin_num = ggplot2::cut_number(!!sym(rt_col), n_bins, labels = FALSE)) %>%
        dplyr::ungroup() -> bins_data

    bins_data %>%
        group_by(dplyr::across(all_of(c(subject_col, factors_col, congruency_col, "bin_num")))) %>%
        dplyr::summarise(mean_RT = mean(!!sym(rt_col), na.rm = TRUE)) %>%
        pivot_wider(
            names_from = !!sym(congruency_col),
            values_from = mean_RT # nolint
        ) %>%
        mutate(simon_effect = !!sym(inc_value) - !!sym(con_value)) %>%
        dplyr::ungroup() -> simon_effect_data

    if (debug) {
        print(simon_effect_data)
    }
    
    # 对Simon effect Data 进行 repeated measures ANOVA
    # anova_no_volatility_results
    anova_no_volatility_results <- rlang::eval_tidy(rlang::expr(ez::ezANOVA(
                data = simon_effect_data,
                dv = .(simon_effect), # 依赖变量
                wid = .(!!sym(subject_col)), # 被试列名
                within = .(prop, bin_num), # 组内变量
                detailed = TRUE # 提供详细的输出
    )))

    # anova_add_volatility_results
    anova_add_volatility_results <- rlang::eval_tidy(rlang::expr(ez::ezANOVA(
                data = simon_effect_data,
                dv = .(simon_effect), # 依赖变量
                wid = .(!!sym(subject_col)), # 被试列名
                within = .(!!!rlang::syms(factors_col), bin_num), # 组内变量
                detailed = TRUE # 提供详细的输出
    )))

    # 对bin_num为1的simon effect数据进行t检验
    simon_effect_data %>%
        filter(bin_num == 1) -> bin1_data
    bin1_MC <- filter(bin1_data, prop == "MC")
    bin1_MI <- filter(bin1_data, prop == "MI")
    ttest_results <- t.test(bin1_MC$simon_effect, bin1_MI$simon_effect, paired = TRUE)


    cat("\033[31m====== The Statistical results of bin1 are as follows: =====\033[39m\n")
    cat("\033[32m============== The t test results ==========================\033[39m\n")
    print(apa::t_apa(ttest_results))
    cat("\033[32m============== The ANOVA results ===========================\033[39m\n")
    print("Without volatility")
    print(anova_no_volatility_results)
    print("Add volatility")
    print(anova_add_volatility_results)

    simon_effect_data %>%
        group_by(dplyr::across(all_of(c(factors_col, "bin_num")))) %>%
        dplyr::summarise(mean_simon_effect = mean(.data$simon_effect,
            na.rm = TRUE
        )) -> delta_plot_data

    ggplot2::ggplot(delta_plot_data, ggplot2::aes(
        x = bin_num,
        y = mean_simon_effect,
        color = !!sym(color_by)
    )) +
        ggplot2::geom_path(size = 1.5) +
        ggplot2::facet_grid(cols = vars(!!sym(facet_by))) +
        ggplot2::ylab("Simon 效应（ms）") +
        ggplot2::xlab("RT 分位数") +
        ggplot2::scale_x_continuous(breaks = 1:n_bins) +
        ggplot2::scale_color_discrete(name = "一致性比例") +
        ggpubr::theme_pubr() +
        ggplot2::theme(
            text = ggplot2::element_text(family = "PingFangSC-Regular", size = 12),
            axis.text = ggplot2::element_text(family = "PingFangSC-Regular", size = 12),
            axis.title = ggplot2::element_text(family = "PingFangSC-Regular", size = 12)
        ) -> delta_plot_pic

    return(list(delta_plot_pic, 
           list(ANOVA_no_volatility = anova_no_volatility_results, 
                ANOVA_with_volatility = anova_add_volatility_results, 
                ttest_bin1 = ttest_results)))
}