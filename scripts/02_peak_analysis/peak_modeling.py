"""
peak_modeling.py - 对peak数据进行分布建模分析

功能：
1. 读取输入的peak数据（feather格式）
2. 数据预处理和过滤
3. 自定义分布模型和正态分布拟合
4. 模型比较（AIC/BIC）
5. 结果保存

使用方法：
python peak_modeling.py --input <输入文件> --output <输出目录> --prefix <输出前缀> [--sample_size 200000]
"""

import argparse
import logging
import os
import json
import numpy as np
import pandas as pd
from scipy.special import erf
from scipy.optimize import minimize
from scipy.stats import norm

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class PeakModeler:
    def __init__(self, config=None):
        """
        初始化模型配置
        
        参数:
        config : dict, 可选
            包含以下键的配置字典:
            - initial_guess: 初始参数猜测 [mu1, sigma1, sigma2, a, b]
            - bounds: 参数优化边界
            - sample_size: 采样数据量
        """
        self.config = config or {}
        self._set_default_config()

    def _set_default_config(self):
        """设置默认配置"""
        default_config = {
            'model_params': {
                'initial_guess': [-100, 90, 10, -500, 0],  
                'bounds': [
                    (None, None),  # mu1
                    (1e-6, None),  # sigma1 > 0
                    (1e-6, None),  # sigma2 > 0
                    (None, -100),  # a
                    (-100, 100),  # b
                ],
                'optimize_method': 'Nelder-Mead',
                'pas_relative_position': 600,
            },
        }
        for k, v in default_config.items():
            self.config.setdefault(k, v)

    def load_data(self, input_path):
        """加载并预处理数据"""
        logger.info(f"Loading data from {input_path}")
        
        # 读取数据
        self.df = pd.read_feather(input_path)
        
        # 数据过滤
        self._preprocess_data()
        
        # 计算fragment length
        self.df["fragment_length"] = self.df["reference_length"] - self.df["distance_to_pas"]
        
        logger.info(f"Total peaks after filtering: {len(self.df):,}")

    def _preprocess_data(self):
        """数据预处理"""
        # 按pas分组过滤
        pas_counts = self.df.groupby("pas")["pas"].transform("count")
        self.df = self.df[pas_counts > 50].reset_index(drop=True)
        num_peaks = self.df["pas"].nunique()
        self.coverage = np.zeros((621, num_peaks), dtype=np.int32)
        for i, (pas, peak_df) in enumerate(self.df.groupby("pas")):
            start_pos = peak_df['distance_to_pas'].values - peak_df['reference_length'].values
            end_pos = peak_df['distance_to_pas'].values
            
            start_index = np.clip(start_pos + 600, 0, 620)
            end_index = np.clip(end_pos + 600, 0, 620)
            
            for start, end in zip(start_index, end_index):
                if start < end:  # 添加有效性检查
                    self.coverage[start:end, i] += 1

    def _likelihood(self, theta, z):
        """自定义似然函数"""
        mu1, sigma1, sigma2, a, b = theta
        sigma_total_sq = sigma1**2 + sigma2**2
        sigma_total = np.sqrt(sigma_total_sq)
        
        numerator = np.exp(-((z - mu1)**2) / (2 * sigma_total_sq))
        
        erf_arg1 = ((a - z) * sigma1**2 + (a - mu1) * sigma2**2) / (
            np.sqrt(2) * sigma1 * sigma2 * sigma_total
        )
        erf_arg2 = ((b - z) * sigma1**2 + (b - mu1) * sigma2**2) / (
            np.sqrt(2) * sigma1 * sigma2 * sigma_total
        )
        erf_diff = erf(erf_arg2) - erf(erf_arg1)
            
        erf_denom1 = (a - mu1) / (np.sqrt(2) * sigma1)
        erf_denom2 = (b - mu1) / (np.sqrt(2) * sigma1)
        erf_denom_diff = erf(erf_denom2) - erf(erf_denom1)
        
        L = (numerator * erf_diff) / (
            np.sqrt(2 * np.pi) * erf_denom_diff * sigma_total
        )
        return np.clip(L, 1e-20, None)

    def _neg_log_likelihood(self, theta, counts, pas_relative_position):
        """
        参数:
        - theta: [mu1, sigma1, sigma2, a, b]
        - counts: 分箱计数数组，索引对应覆盖度值
        - pas_relative_position: pas 的相对位置
        """
        # 计算所有分箱的概率密度
        z_values = np.arange(len(counts)) - pas_relative_position
        pdfs = self._likelihood(theta, z_values)
        
        # 转换为概率（因分箱宽度=1，概率=密度×1）
        probs = pdfs
        
        # 计算对数似然
        valid_bins = counts > 0
        log_lik = np.sum(counts[valid_bins] * np.log(probs[valid_bins]))
        return -log_lik
    
    def weighted_normal_negloglik(self, params, bin_centers, counts):
        mu, sigma = params
        # 计算每个分箱的概率（精确积分）
        probs = norm.cdf(bin_centers + 0.5, mu, sigma) - norm.cdf(bin_centers - 0.5, mu, sigma)
        probs = np.clip(probs, 1e-20, 1.0)
        return -np.sum(counts * np.log(probs))

    def fit_models(self):
        """拟合所有模型"""
        # 采样数据
        self.z_data = self.coverage.sum(axis=1)
        # 拟合自定义模型
        logger.info("Fitting custom model...")
        result = minimize(
            self._neg_log_likelihood,
            self.config['model_params']['initial_guess'],
            args=(self.z_data, self.config['model_params']['pas_relative_position']),
            method=self.config['model_params']['optimize_method'],
            bounds=self.config['model_params']['bounds'],
            options={
                "maxiter": 5000
            }
        )
        
        if not result.success:
            logger.error("Custom model fitting failed!")
            raise RuntimeError(result.message)
            
        self.custom_params = result.x
        self.custom_loglik = -result.fun
        # self.custom_loglik = -self._neg_log_likelihood(self.custom_params, self.df["fragment_length"].values)
        
        # 拟合正态分布

        # 初始猜测值
        bin_centers = np.arange(len(self.z_data)) - self.config['model_params']['pas_relative_position']
        initial_mu = np.average(bin_centers, weights=self.z_data)
        initial_sigma = np.sqrt(np.average((bin_centers - initial_mu)**2, weights=self.z_data))

        # 优化拟合
        result_norm = minimize(
            self.weighted_normal_negloglik,
            [initial_mu, initial_sigma],
            args=(bin_centers, self.z_data),
            bounds=[(None, None), (1e-6, None)]
        )
        mu_norm, sigma_norm = result_norm.x

        logger.info("Fitting normal distribution...")
        self.norm_params = result_norm.x
        self.norm_loglik = -result_norm.fun
        # self.norm_loglik = np.sum(norm.logpdf(self.df["fragment_length"].values, *self.norm_params))
        
        return self

    def evaluate_models(self):
        """模型评估"""
        n = int(self.coverage.sum())
        
        # 自定义模型（参数数量5）
        aic_custom = 2*5 - 2*self.custom_loglik
        bic_custom = 5*np.log(n) - 2*self.custom_loglik
        
        # 正态分布（参数数量2）
        aic_norm = 2*2 - 2*self.norm_loglik
        bic_norm = 2*np.log(n) - 2*self.norm_loglik
        
        apex_positions = []
        weighted_positions = []
        for i in range(self.coverage.shape[1]):
            # 获取当前peak的覆盖数据列
            column = self.coverage[:, i]
            
            # 找到最大值及其所有位置索引
            max_val = np.max(column)
            max_indices = np.where(column == max_val)[0]
            
            # 计算平均索引并转换为相对于PAS的位置
            avg_index = np.mean(max_indices)
            apex_relative = avg_index - 600  # 根据数据结构调整偏移量
            
            apex_positions.append(apex_relative)

            column_indices = np.arange(len(column))
            weighted_position = np.sum(column * column_indices) / np.sum(column) - 600
            weighted_positions.append(weighted_position)
            
        apex_positions = np.array(apex_positions)
        weighted_positions = np.array(weighted_positions)
        apex_mean = float(np.mean(apex_positions))
        weighted_mean = float(np.mean(weighted_positions))
        apex_std = float(np.std(apex_positions, ddof=1))
        weighted_std = float(np.std(weighted_positions, ddof=1))

        self.results = {
            "custom_model": {
                "params": {
                    "mu1": self.custom_params[0],
                    "sigma1": self.custom_params[1],
                    "sigma2": self.custom_params[2],
                    "a": self.custom_params[3],
                    "b": self.custom_params[4]
                },
                "log_likelihood": self.custom_loglik,
                "AIC": aic_custom,
                "BIC": bic_custom
            },
            "normal_dist": {
                "params": {
                    "mu": self.norm_params[0],
                    "sigma": self.norm_params[1]
                },
                "log_likelihood": self.norm_loglik,
                "AIC": aic_norm,
                "BIC": bic_norm
            },
            "data_size": n,
            "apex_mean": apex_mean,
            "weighted_mean": weighted_mean,
            "apex_std": apex_std,
            "weighted_std": weighted_std
        }
        return self

    def save_results(self, output_dir, prefix):
        """保存结果"""
        os.makedirs(output_dir, exist_ok=True)
        
        # 保存参数和指标
        with open(os.path.join(output_dir, f"{prefix}_model_results.json"), "w") as f:
            json.dump(self.results, f, indent=2)
        
        # # 保存采样数据
        # np.save(os.path.join(output_dir, "sampled_data.npy"), self.z_data)
        
        logger.info(f"Results saved to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='peak片段长度分布建模')
    parser.add_argument('--input', required=True, help='输入peak数据路径（feather格式）')
    parser.add_argument('--output', required=True, help='输出目录')
    parser.add_argument('--prefix', required=True, help='输出文件前缀')
    args = parser.parse_args()

    try:
        # 初始化建模器
        modeler = PeakModeler()
        
        # 加载和处理数据
        modeler.load_data(args.input)
        
        # 模型拟合
        modeler.fit_models()
        
        # 模型评估
        modeler.evaluate_models()
        
        # 保存结果
        print(modeler.results)
        modeler.save_results(args.output, args.prefix)
        
        # 打印结果摘要
        print("\n模型比较结果:")
        print(f"{'模型':<15} {'AIC':<15} {'BIC':<15}")
        print(f"{'自定义模型':<15} {modeler.results['custom_model']['AIC']:<15.2f} {modeler.results['custom_model']['BIC']:<15.2f}")
        print(f"{'正态分布':<15} {modeler.results['normal_dist']['AIC']:<15.2f} {modeler.results['normal_dist']['BIC']:<15.2f}")
        
    except Exception as e:
        logger.error(f"处理失败: {str(e)}")
        raise

if __name__ == "__main__":
    main()
