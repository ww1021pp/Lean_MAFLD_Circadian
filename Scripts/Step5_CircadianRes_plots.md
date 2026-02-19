rm(list=ls())
setwd("~/project/ProtMetabolite_PredictCT/UKBiobank_data/MAFLD_definition/Circadian_diff_LeanMAFLD/HarmReg_res/")

library(dplyr)
library(data.table)
library(tidyr)

res.norm.df.list = readRDS("NormE.HarReg_resForLeanMAFLD_downSample.rds")

fit_para_controls_norm <- bind_rows(lapply(res.norm.df.list[1:10], function(x){
  df <- x$fit.para
  df$group <- "Control"
  df
}))

# Add group to case
fit_para_case_norm <- res.norm.df.list$MAFLD$fit.para
fit_para_case_norm$group <- "Case"

#####Aggregate across runs for stability
summary_ctrl_norm <- fit_para_controls_norm %>%
  filter(group == "Control") %>%
  group_by(Gene) %>%
  summarise(
    mean_log2_pt = mean(amp.c),
    sd_log2_pt = sd(amp.c),
    mean_phase = mean(phase.s),
    sd_phase = sd(phase.s),
    max_BHQ = max(BHQ, na.rm = TRUE), 
    .groups = "drop"
  )

##### Define “stable” circadian genes
####Set thresholds for stability:####

stable_genes_norm <- summary_ctrl_norm %>%
  filter(
    sd_log2_pt < 0.1,   # amplitude is consistent
    sd_phase < 2,        # phase is consistent (hours)
    mean_log2_pt > 0.25,   # biologically meaningful amplitude
    )
library(dplyr)
library(circular)
library(ggplot2)

# --- 1. PREPARE DATASETS ---
# Ensure your columns are named consistently
# Ctrl: protein, mean_log2_pt (amp), mean_phase
# Case: protein, amp.c (amp), phase_case, BHQ

ctrl_all <- summary_ctrl_norm %>%
  select(Gene, amp_ctrl = mean_log2_pt, phase_ctrl = mean_phase,BHQ_ctrl =max_BHQ)

case_all <- fit_para_case_norm %>%
  select(Gene, amp_case = amp.c, phase_case = phase.s, BHQ_case = BHQ)

Ctrl_cir = stable_genes_norm %>%
           select(Gene, amp_ctrl = mean_log2_pt, phase_ctrl = mean_phase,BHQ_ctrl =max_BHQ)

Case_cir = fit_para_case_norm %>% filter(BHQ <1e-5, amp.c>0.25) %>%
  select(Gene, amp_case = amp.c, phase_case = phase.s, BHQ_case = BHQ)


# --- 2. JOIN AND CATEGORIZE ---
all_background <- full_join(ctrl_all, case_all, by = "Gene") %>%
  mutate(category = case_when(
    (Gene %in% Ctrl_cir$Gene) & (Gene %in% Case_cir$Gene)  ~ "Common",
    (Gene %in% Ctrl_cir$Gene) & !(Gene %in% Case_cir$Gene)   ~ "Ctrl-Specific",
    !(Gene %in% Ctrl_cir$Gene) & (Gene %in% Case_cir$Gene)   ~ "Case-Specific",
    TRUE ~ "non_circadian"
  ))


# --- 3. PHASE SHIFT ANALYSIS (For Common Group Only) ---
# 2. Join and Categorize using your defined Filtered Lists
full_comparison <- all_background %>%
  mutate(category = case_when(
    (Gene %in% Ctrl_cir$Gene) & (Gene %in% Case_cir$Gene)  ~ "Common",
    (Gene %in% Ctrl_cir$Gene) & !(Gene %in% Case_cir$Gene) ~ "Ctrl-Specific",
    !(Gene %in% Ctrl_cir$Gene) & (Gene %in% Case_cir$Gene) ~ "Case-Specific",
    TRUE ~ "non_circadian"  # Corrected syntax for 'else'
  ))

# 3. PHASE SHIFT ANALYSIS (Only for Common Proteins)
# We calculate the shortest distance on a 24h circle
# 3. PHASE SHIFT ANALYSIS (Only for Common Proteins)
# We calculate the shortest distance on a 24h circle
common_analysis <- full_comparison %>%
  filter(category == "Common") %>%
  mutate(
    # Convert to radians for circular math
    rad_ctrl = ( phase_ctrl / 24) * 2 * pi,
    rad_case = (phase_case / 24) * 2 * pi,
    
    # Circular difference: atan2 handles the 23:00 to 01:00 jump correctly
    diff_rad = atan2(sin(rad_case - rad_ctrl), cos(rad_case - rad_ctrl)),
    
    # Convert back to hours (-12 to 12)
    phase_shift_hrs = diff_rad * (12 / pi)
  )


common_diff_analysis <- common_analysis %>%
  mutate(
    # 1. Amplitude Change (Log2 Fold Change)
    amp_l2fc = log2(amp_case / amp_ctrl),
    
    # 2. Categorize by both Phase and Amp
    # We use 1 hour as a significant phase shift and 0.25 for L2FC
    diff_category = case_when(
      amp_l2fc > 0.25 & phase_shift_hrs > 1   ~ "Stronger & Delayed",
      amp_l2fc > 0.25 & phase_shift_hrs < -1  ~ "Stronger & Advanced",
      amp_l2fc < -0.25 & phase_shift_hrs > 1  ~ "Weaker & Delayed",
      amp_l2fc < -0.25 & phase_shift_hrs < -1 ~ "Weaker & Advanced",
      abs(amp_l2fc) < 0.25 & abs(phase_shift_hrs) < 1 ~ "Stable",
      TRUE ~ "Mixed Change"
    )
  )
fwrite(x = common_diff_analysis,"Common_circadian_genes_diffRadAndAmp.txt",sep = "\t",quote = F,row.names = F)


############### plot phase distribution########3
# 1. Create the 'Control' Population (Common + Ctrl-Specific)
control_population <- full_comparison %>%
  filter(category %in% c("Common", "Ctrl-Specific")) %>%
  select(Gene, phase = phase_ctrl,amp=amp_ctrl,BHQ=BHQ_ctrl,category=category) %>%
  mutate(Group = "Health")

# 2. Create the 'Case' Population (Common + Case-Specific)
case_population <- full_comparison %>%
  filter(category %in% c("Common", "Case-Specific")) %>%
  select(Gene, phase = phase_case,amp=amp_case,BHQ=BHQ_case,category=category) %>%
  mutate(Group = "LeanMAFLD")
fwrite(control_population,"All_control_Circadian_genes.txt",quote = F,sep="\t",row.names = F)
fwrite(control_population,"All_case_Circadian_genes.txt",quote = F,sep="\t",row.names = F)


# 3. Combine them
comparison_bar_data <- bind_rows(control_population, case_population)

###### plot phase distribution#####
ggplot(comparison_bar_data,aes(x = phase, fill = Group)) +
  geom_histogram(binwidth = 0.8,colour = "black",position = "Identity",alpha=0.6) +
  scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 24)) +
  scale_fill_manual(values = c("Health" = "#377eb8", "LeanMAFLD" = "#e41a1c")) +
  # --- THE STYLING ---
  theme_bw() + # Start with Black & White theme
  theme(
    panel.grid.major = element_blank(), # Remove major grid
    panel.grid.minor = element_blank(), # Remove minor grid
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1), # White bg, Black box
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1), # Ensure the box is solid
    axis.line = element_blank(), # Let the box serve as the axes
    legend.background = element_blank()
  ) +
  labs(
    x = "Phase(h)",
    y = "Protein Count"
  )
ggsave("../Figs_out/Distribution_of_phase.pdf",height = 3,width = 4)
<img width="458" height="347" alt="image" src="https://github.com/user-attachments/assets/2ac310aa-7afe-42b2-a0ca-8f29659b2d68" />




###################### plot common , Control Specific and Case specific phase distribution rose plot ###

make_rose <- function(data, cat_name, fill_color) {
  ggplot(subset(data, category == cat_name), aes(x = phase)) +
    geom_histogram(breaks = seq(0, 24, by = 2),fill = fill_color,colour = "black") +
    coord_polar(start = 0) +
    scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24), expand = c(0, 0)) +
    theme_minimal() +
    labs(x = cat_name, y = "Count") +
    theme(text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 8,colour = "black"),
          axis.text = element_text(size = 8,colour = "black"))

}

# Create the individual plots
p1 <- make_rose(comparison_bar_data, "Common", "grey60")
p2 <- make_rose(comparison_bar_data, "Ctrl-Specific", "#377EB8")
p3 <- make_rose(comparison_bar_data, "Case-Specific", "#E41A1C")

# Combine them in a row (1 row, 3 columns)
three_row_plot <- p2 + p1 + p3 + plot_layout(nrow = 1)
ggsave("../Figs_out/RosePlot_of_phasefor3Catogory.pdf",height = 3,width = 6)

<img width="702" height="264" alt="image" src="https://github.com/user-attachments/assets/85ce2b69-433b-49fe-bf60-553580b7dbac" />


###### plot 
make_density <- function(data, cat_name, fill_color) {
  ggplot(subset(data, category == cat_name), aes(x = phase))  +
   geom_density(fill = fill_color,  color = "black", linewidth = 0.5) +
       scale_x_continuous(breaks = seq(0, 24, 2), limits = c(0, 24)) +
       labs(title = cat_name ,
            x = "Phase (ZT)", y = "Density of Proteins")+
       theme_bw() +
       theme(text = element_text(size = 8, colour = "black"),
          axis.title = element_text(size = 8,colour = "black"),
          axis.text = element_text(size = 8,colour = "black"),
          panel.grid.major = element_blank(), # Remove major grid
          panel.grid.minor = element_blank(), # Remove minor grid
          panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5), # White bg, Black box
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), # Ensure the box is solid
          axis.line = element_blank(), # Let the box serve as the axes
          legend.background = element_blank())
}

p1_1 <- make_density(comparison_bar_data, "Common", "grey60")
p2_1 <- make_density(comparison_bar_data, "Ctrl-Specific", "#377EB8")
p3_1 <- make_density(comparison_bar_data, "Case-Specific", "#E41A1C")


p_density = p2_1 + p1_1 + p3_1 + plot_layout(nrow = 1)

ggsave("../Figs_out/densityDistribution_of_phase_3Catogory.pdf",plot = p_density,height = 2,width = 6)

<img width="700" height="221" alt="image" src="https://github.com/user-attachments/assets/89ba9db3-7e95-4451-8365-87db0aa5e232" />


########################## plot circadian phase amplitude scatter plot###
Top_hits = control_population %>% filter(category=="Ctrl-Specific") %>% arrange(BHQ) %>% slice_head(n=30)

library(ggrepel)
p1_ctrl = ggplot(control_population, aes(x = phase, y = amp, color = category)) +
  geom_point(     data = subset(control_population, category == "Common"),
                  color = "grey60",size = 0.1) +
  geom_point(data = subset(control_population, category == "Ctrl-Specific"),
             aes(color = -log10(BHQ)),size = 0.3)+
  geom_text(
    data = Top_hits, 
    aes(label = Gene), colour = "black",
    size = 3, check_overlap = TRUE,
    vjust = -1 # Puts the text slightly above the dot
  ) +scale_color_viridis_c(option = "C", name = "-log10(BHQ)")+
  scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 24)) +
  theme_bw() +
 # theme(panel.grid = element_blank())+ # Start with Black & White theme
  theme(
    panel.grid.major = element_blank(), # Remove major grid
    panel.grid.minor = element_blank(), # Remove minor grid
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1), # White bg, Black box
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1), # Ensure the box is solid
    axis.line = element_blank(), # Let the box serve as the axes
    legend.background = element_blank()
  ) +
  labs(
    x = "Phase(h)",
    y = "amplitude(log2(FC))"
  )

Top_case = case_population %>% filter(category == "Case-Specific") %>% arrange(BHQ) %>% head(30)
p1_case = ggplot(case_population, aes(x = phase, y = amp, color = category)) +
  geom_point(     data = subset(case_population, category == "Common"),
                  color = "grey60",size = 0.1) +
  geom_point(data = subset(case_population, category == "Case-Specific"),
             aes(color = -log10(BHQ)),size = 0.3)+
  geom_text(
    data = Top_case, 
    aes(label = Gene), colour = "black",
    size = 3, check_overlap = TRUE,
    vjust = -1 # Puts the text slightly above the dot
  ) +scale_color_viridis_c(option = "C", name = "-log10(BHQ)")+
  scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 24)) +
  theme_bw() +
  # theme(panel.grid = element_blank())+ # Start with Black & White theme
  theme(
    panel.grid.major = element_blank(), # Remove major grid
    panel.grid.minor = element_blank(), # Remove minor grid
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1), # White bg, Black box
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1), # Ensure the box is solid
    axis.line = element_blank(), # Let the box serve as the axes
    legend.background = element_blank()
  ) +
  labs(
    x = "Phase(h)",
    y = "amplitude(log2(FC))"
  )

library(patchwork)
p= p1_ctrl + p1_case
ggsave(filename = "../Figs_out/scatterplot_phase_amp_control_MALFD.pdf",plot = p,height = 4, width = 10)

<img width="1174" height="458" alt="image" src="https://github.com/user-attachments/assets/2a54778f-6e02-40f8-a283-2c9db56f3941" />


##############################
library(circular)
library(ggplot2)

# 1. Calculate the Circular Correlation Coefficient
# Uses the 'Common' group defined in the previous step
cor_result <- cor.circular(
  circular(common_analysis$rad_ctrl, units="radians"),
  circular(common_analysis$rad_case, units="radians")
)

print(paste("Circular Correlation (rho):", round(cor_result, 3)))

# 2. Visualize the Phase Alignment
ggplot(common_analysis, aes(x = phase_ctrl, y = phase_case)) +
  geom_point(alpha = 0.5, color = "midnightblue") +
  # Add a diagonal line representing 'Perfect Alignment'
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  # Add a line representing the 'Systemic Shift' if applicable
  #geom_abline(slope = 1, intercept = as.numeric(global_shift), 
  #            linetype = "dotted", color = "blue") +
  theme_bw() + # Start with Black & White theme
  theme(
    panel.grid.major = element_blank(), # Remove major grid
    panel.grid.minor = element_blank(), # Remove minor grid
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1), # White bg, Black box
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1), # Ensure the box is solid
    axis.line = element_blank(), # Let the box serve as the axes
    legend.background = element_blank()
  ) +
  coord_fixed(xlim = c(0, 24), ylim = c(0, 24)) +
  labs(
    title = "Phase Alignment: Control vs. Case",
    subtitle = paste("Circular Correlation:", round(cor_result, 3)),
    x = "Control Peak Time (h)",
    y = "Case Peak Time (h)"
  )
ggsave("../Figs_out/Correlation_plot_of_Scatter_plots_CtrlphaseVSCasephase.pdf",height =4,width = 4 )
<img width="461" height="466" alt="image" src="https://github.com/user-attachments/assets/95ce404f-dc8c-4e5d-a139-a4f17d5a0ecf" />

