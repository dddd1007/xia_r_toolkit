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
        ggplot2::xlab("反应时分位数") +
        ggplot2::scale_x_continuous(breaks = 1:n_bins) +
        ggplot2::scale_color_discrete(name = "一致性比例") +
        ggpubr::theme_pubr() +
        ggplot2::theme(
            text = ggplot2::element_text(family = "PingFangSC-Regular", size = 12),
            axis.text = ggplot2::element_text(family = "PingFangSC-Regular", size = 12),
            axis.title = ggplot2::element_text(family = "PingFangSC-Regular", size = 12)
        ) -> delta_plot_pic

    return(list(
        delta_plot = delta_plot_pic,
        simon_effect_data = simon_effect_data,
        delta_plot_data = delta_plot_data
    ))
}

#' rt_dist_anova_and_ttest
#'
#' 计算重复测量ANOVA和事后检验
#'
#' @param simon_effect_data Simon效应数据，应为一个数据框
#' @param subject_col 被试列名
#' @param factors_col 因子列名，应为一个字符向量
#' @param print_result 是否打印结果，默认为 TRUE
#' @return 返回一个列表，包含ANOVA结果，描述统计结果和t检验结果
#' @export
#' @examples
#' # 假设我们有一个名为 simon_effect_data 的数据框，其中包含因子 "subject_num" 和 "condition"
#' rt_dist_anova_and_ttest(simon_effect_data,
#'     subject_col = "subject_num",
#'     factors_col = c("condition"), print_result = TRUE
#' )
rt_dist_anova_and_ttest <- function(simon_effect_data,
                                    subject_col,
                                    factors_col,
                                    facet_col,
                                    print_result = TRUE) {
    simon_effect_data %>%
        group_by(!!!syms(factors_col), bin_num) %>%
        get_summary_stats(simon_effect, type = "mean_sd") -> all_factor_summary
    all_factor_summary %>%
        select(-n, -sd) %>%
        pivot_wider(names_from = prop, values_from = mean) %>%
        mutate(prop_delta = abs(MI - MC)) %>%
        ungroup() -> prop_delta
    prop_delta %>%
        select(-MC, -MI) %>%
        group_by(bin_num) %>%
        summarise(mean_prop_delta = mean(prop_delta)) %>%
        ungroup() -> mean_prop_delta

    # 使用rstatix包进行重复测量ANOVA
    anova_add_volatility_results <- simon_effect_data %>%
        anova_test(
            dv = simon_effect,
            wid = !!sym(subject_col),
            within = c(all_of(factors_col), "bin_num")
        )

    # 进行事后检验
    # Simple two-way interaction
    simon_effect_data %>%
        group_by(bin_num, !!sym(facet_col)) %>%
        anova_test(dv = simon_effect, wid = !!sym(subject_col), within = prop) -> between_prop
    simon_effect_data %>%
        group_by(!!sym(facet_col)) %>%
        anova_test(dv = simon_effect, wid = !!sym(subject_col), within = c(bin_num, prop)) -> between_prop_bin
    simon_effect_data %>%
        group_by(prop, !!sym(facet_col)) %>%
        anova_test(
            dv = simon_effect,
            wid = !!sym(subject_col),
            within = bin_num
        ) %>%
        get_anova_table() -> between_bin

    # simple simple effects
    posthoc_results <- simon_effect_data %>%
        group_by(!!!syms(factors_col)) %>%
        pairwise_t_test(simon_effect ~ bin_num, paired = TRUE, p.adjust.method = "bonferroni")

    # 对bin_num为1的simon effect数据进行t检验
    simon_effect_data %>%
        filter(bin_num == 1) -> bin1_data
    bin1_MC <- filter(bin1_data, prop == "MC")
    bin1_MI <- filter(bin1_data, prop == "MI")
    ttest_results <- t.test(bin1_MC$simon_effect, bin1_MI$simon_effect, paired = TRUE)

    # 输出描述统计结果
    if (print_result) {
        line_length <- 50 # 您可以根据需要调整这个数字
        # 计算中文字符数目
        count_chinese <- function(string) {
            n <- nchar(string, type = "bytes") - nchar(string, type = "chars")
            return(n / 2)
        }

        # 格式化输出行
        format_line <- function(text, total_length, color_code = "32") {
            chinese_count <- count_chinese(text)
            side_length <- (total_length - nchar(text) - 2 - chinese_count) %/% 2
            line <- paste0(
                paste0("\033[", color_code, "m"),
                strrep("=", side_length), " ", text, " ",
                strrep("=", side_length), "\033[39m\n"
            )
            return(line)
        }

        cat(format_line("开始输出统计分析报告", line_length, "31"))
        cat(format_line("描述统计结果", line_length))
        cat("\033[34mAll factor summary:\033[39m\n")
        print(all_factor_summary)
        cat("\033[34mBetween prop delta:\033[39m\n")
        print(prop_delta)
        cat("\033[34mMean prop delta:\033[39m\n")
        print(mean_prop_delta)

        cat(format_line("ANOVA结果", line_length))
        print(anova_add_volatility_results)

        # 输出事后检验结果
        cat(format_line("事后检验结果", line_length))
        cat("\033[34mSimple two-way interaction:\033[39m\n")
        print(between_prop)
        cat("\033[34mSimple Simple Effect:\033[39m\n")
        print(posthoc_results)

        # 输出统计结果
        cat(format_line("Bin 1 ttest 统计结果", line_length))
        print(report::report(ttest_results))

        cat(format_line("统计报告分析结果输出结束", line_length, "31"))
    }


    return(list(
        anova_results = anova_add_volatility_results,
        desc_stats = list(
            all_factor_summary = all_factor_summary,
            prop_delta = prop_delta, mean_prop_delta = mean_prop_delta
        ),
        ttest_results = ttest_results
    ))
}
