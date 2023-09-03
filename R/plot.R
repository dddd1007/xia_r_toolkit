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
                          print_desc = TRUE) {
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

    if (print_desc) {
        # 将 tibble data_forplot 转换为 dataframe
        data_forplot <- as.data.frame(data_forplot)
        print(data_forplot)
        input_data %>%
            as.data.frame() %>%
            dplyr::group_by(!!!rlang::syms(setdiff(factors, facet_by))) %>%
            dplyr::summarise(mean_mean_dep_var = sprintf("%.4f", mean(!!sym(dep_var)))) %>%
            print()
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
