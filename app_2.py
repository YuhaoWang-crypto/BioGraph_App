import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import re
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
# 这里的 CountVectorizer 如果没用到可以去掉，减小内存占用

# ==========================================
# 1. 页面配置
# ==========================================
st.set_page_config(
    page_title="BioGraph: Protein Lifecycle Explorer",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🧬 BioGraph: 多维蛋白质组学流形分析平台")
st.markdown("""
**真实数据驱动**：本平台展示了基于 UniProt 真实序列与模拟组学数据的深度分析结果。
已预先集成了 UMAP 流形投影、K-Means 聚类及功能富集分析。
""")

# ==========================================
# 2. 数据加载 (读取 GZ 压缩文件)
# ==========================================
@st.cache_data
def load_data():
    """
    读取压缩的 CSV 文件 (.csv.gz)。
    Pandas 会根据文件后缀自动推断解压方式，或者显式指定 compression='gzip'。
    """
    try:
        # === 修改点：读取 .csv.gz 文件 ===
        # Pandas 能够自动识别 .gz 后缀并解压
        df = pd.read_csv("final_analysis_result.csv.gz", compression='gzip')
        
        # --- 下面是常规清洗逻辑 ---
        
        # 简单的清洗，防止 NaN 报错
        if 'cc_function' in df.columns:
            df['cc_function'] = df['cc_function'].fillna('Unknown')
        
        if 'Gene_Symbol' in df.columns:
            df['Gene_Symbol'] = df['Gene_Symbol'].fillna('Unknown')
        
        # 确保半衰期数值正常
        if 'Real_Protein_HalfLife_Hours' in df.columns:
            df['Real_Protein_HalfLife_Hours'] = df['Real_Protein_HalfLife_Hours'].fillna(0)
            
        # 生成辅助标签
        if 'Real_Protein_HalfLife_Hours' in df.columns:
            df['Stability_Level'] = pd.cut(df['Real_Protein_HalfLife_Hours'], 
                                           bins=[-1, 10, 50, 10000], 
                                           labels=['Short (<10h)', 'Medium', 'Long (>50h)']).astype(str)
        
        return df
    except FileNotFoundError:
        st.error("❌ 未找到数据文件！请确保 'final_analysis_result.csv.gz' 已上传到项目根目录。")
        return pd.DataFrame()
    except Exception as e:
        st.error(f"❌ 数据读取出错: {e}")
        return pd.DataFrame()

# 加载主数据
df_main = load_data()

if df_main.empty:
    st.stop() # 如果没数据，停止运行

# ==========================================
# 3. 快速重算 PCA Loadings (为了画箭头图)
# ==========================================
@st.cache_data
def calculate_pca_loadings(df):
    """
    因为 CSV 里只有坐标没有模型，这里快速重跑一次 PCA 以获取因子载荷。
    """
    # 选取数值列 (必须与 Colab 分析时一致)
    features = ['Real_Protein_HalfLife_Hours', 'mRNA_Expression', 'Circadian_Amplitude']
    
    # 确保列存在
    valid_features = [f for f in features if f in df.columns]
    
    if len(valid_features) < 2:
        return None, None
        
    # 取对数并标准化
    X = np.log1p(df[valid_features].fillna(0))
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=2)
    pca.fit(X_scaled)
    
    loadings = pd.DataFrame(
        pca.components_.T, 
        columns=['PCA1', 'PCA2'], 
        index=valid_features
    )
    return loadings, valid_features

df_loadings, used_features = calculate_pca_loadings(df_main)

# ==========================================
# 4. 侧边栏：全局过滤器
# ==========================================
st.sidebar.header("🔍 数据筛选")

# 1. Cluster 筛选
if 'Cluster' in df_main.columns:
    all_clusters = sorted(df_main['Cluster'].unique())
    selected_clusters = st.sidebar.multiselect(
        "筛选 Cluster ID",
        options=all_clusters,
        default=all_clusters
    )
    # 应用筛选
    df_filtered = df_main[df_main['Cluster'].isin(selected_clusters)]
else:
    df_filtered = df_main

# 2. 搜索基因高亮
search_gene = st.sidebar.text_input("搜索特定基因 (如 TP53, ALB)", "").upper()

st.sidebar.markdown("---")
st.sidebar.info(f"当前展示: {len(df_filtered)} / {len(df_main)} 个蛋白")


# ==========================================
# 5. 主界面 Tabs
# ==========================================
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🌌 全景流形 (UMAP)", 
    "🔍 蛋白雷达 (Search)", 
    "🕸️ 功能网络 (Network)", 
    "📉 PCA 解密 (Loadings)",
    "🧪 专业富集 (Heatmap)"
])

# --- Tab 1: UMAP 全景 ---
with tab1:
    col1, col2 = st.columns([1, 4])
    with col1:
        # 动态获取可用的上色选项
        available_cols = [c for c in ['Cluster', 'Stability_Level', 'N_Term_AA', 'Is_Cancer'] if c in df_filtered.columns]
        # 如果 Is_Cancer 还没计算，且有 cc_function，则算一下
        if 'Is_Cancer' not in df_filtered.columns and 'cc_function' in df_filtered.columns:
             df_filtered['Is_Cancer'] = df_filtered['cc_function'].str.contains('cancer|tumor', case=False).map({True:'Yes', False:'No'})
             available_cols.append('Is_Cancer')
             
        color_col = st.radio("上色依据:", available_cols, index=0)

    with col2:
        if 'UMAP_X' in df_filtered.columns:
            fig = px.scatter(
                df_filtered, 
                x='UMAP_X', y='UMAP_Y', 
                color=color_col,
                hover_data=['Gene_Symbol', 'Real_Protein_HalfLife_Hours', 'cc_function'],
                title=f"Functional Manifold (Colored by {color_col})",
                height=650,
                template="plotly_white",
                opacity=0.6,
                color_discrete_sequence=px.colors.qualitative.Bold
            )
            
            # 标记搜索的基因
            if search_gene and not df_filtered[df_filtered['Gene_Symbol'] == search_gene].empty:
                row = df_filtered[df_filtered['Gene_Symbol'] == search_gene].iloc[0]
                fig.add_trace(go.Scatter(
                    x=[row['UMAP_X']], y=[row['UMAP_Y']],
                    mode='markers+text',
                    marker=dict(size=20, color='red', symbol='star', line=dict(width=2, color='black')),
                    text=[search_gene],
                    textposition="top center",
                    name='Searched'
                ))

            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("数据中缺少 UMAP 坐标列，无法绘图。")

# --- Tab 2: 详情雷达 ---
with tab2:
    gene_list = sorted(df_main['Gene_Symbol'].unique())
    default_idx = 0
    if search_gene and search_gene in gene_list:
        default_idx = gene_list.index(search_gene)
        
    selected_gene_tab2 = st.selectbox("选择基因查看详情:", gene_list, index=default_idx)
    
    if selected_gene_tab2:
        row = df_main[df_main['Gene_Symbol'] == selected_gene_tab2].iloc[0]
        
        c1, c2 = st.columns([1, 1])
        with c1:
            st.subheader(f"🧬 {selected_gene_tab2}")
            st.write(f"**N端氨基酸:** {row.get('N_Term_AA', 'N/A')}")
            st.write(f"**处理类型:** {row.get('Processing_Type', 'N/A')}")
            st.metric("真实半衰期", f"{row.get('Real_Protein_HalfLife_Hours', 0):.1f} h")
            st.metric("mRNA 表达量", f"{row.get('mRNA_Expression', 0):.2f}")
            st.info(f"所属 Cluster: {row.get('Cluster', 'N/A')}")
            
        with c2:
            st.markdown("### 功能描述")
            st.caption(row.get('cc_function', 'No description available.'))
            
            st.markdown("### 坐标信息")
            st.json({
                "UMAP": [round(row.get('UMAP_X', 0), 2), round(row.get('UMAP_Y', 0), 2)],
                "PCA": [round(row.get('PCA1', 0), 2), round(row.get('PCA2', 0), 2)]
            })

# --- Tab 3: 网络 (简化版) ---
with tab3:
    st.markdown("### 动态构建功能模块网络")
    keyword = st.text_input("输入关键词:", "mitochondrial")
    
    if keyword:
        subset = df_main[df_main['cc_function'].str.contains(keyword, case=False, na=False)].head(150)
        
        if len(subset) > 1:
            G = nx.Graph()
            genes = subset['Gene_Symbol'].tolist()
            hls = subset['Real_Protein_HalfLife_Hours'].tolist()
            
            for i in range(len(genes)):
                G.add_node(genes[i])
                for j in range(i+1, len(genes)):
                    if abs(hls[i] - hls[j]) < 2.0:
                        G.add_edge(genes[i], genes[j])
            
            fig_net, ax = plt.subplots(figsize=(10, 8))
            pos = nx.spring_layout(G, k=0.15, seed=42)
            nx.draw_networkx_nodes(G, pos, node_size=30, node_color='teal', alpha=0.6, ax=ax)
            nx.draw_networkx_edges(G, pos, alpha=0.1, ax=ax)
            if len(genes) < 50: 
                nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)
            
            ax.set_title(f"Network for '{keyword}' ({len(genes)} nodes)")
            ax.axis('off')
            st.pyplot(fig_net)
        else:
            st.warning("找不到足够的相关蛋白，请更换关键词。")

# --- Tab 4: PCA Loadings ---
with tab4:
    if df_loadings is not None:
        col_l, col_r = st.columns(2)
        with col_l:
            st.write("### 因子载荷表")
            st.dataframe(df_loadings.style.background_gradient(cmap='RdBu'))
        with col_r:
            st.write("### 向量贡献图")
            fig_l, ax = plt.subplots(figsize=(6, 6))
            ax.axhline(0, color='grey', linestyle='--')
            ax.axvline(0, color='grey', linestyle='--')
            for i, feat in enumerate(df_loadings.index):
                ax.arrow(0, 0, df_loadings.iloc[i, 0], df_loadings.iloc[i, 1], 
                         color='red', width=0.01, head_width=0.05)
                ax.text(df_loadings.iloc[i, 0]*1.2, df_loadings.iloc[i, 1]*1.2, feat, color='darkblue')
            ax.set_xlabel("PCA 1")
            ax.set_ylabel("PCA 2")
            st.pyplot(fig_l)
    else:
        st.warning("无法计算 PCA，可能缺少数值列。")

# --- Tab 5: 富集热图 ---
with tab5:
    st.write("### 专业生物学术语富集度")
    
    BIO_DICT = {
        'Loc': ['mitochondrion', 'nucleus', 'membrane', 'cytoplasm', 'secreted'],
        'Func': ['kinase', 'transcription', 'transport', 'metabolism', 'immune'],
        'Struct': ['zinc', 'finger', 'domain', 'repeat']
    }
    
    keywords = [k for v in BIO_DICT.values() for k in v]
    
    if 'Cluster' in df_main.columns:
        clusters = sorted(df_main['Cluster'].unique())
        heatmap_data = []
        for k in keywords:
            row_data = []
            for c in clusters:
                sub = df_main[df_main['Cluster'] == c]
                if len(sub) > 0:
                    ratio = sub['cc_function'].str.contains(k, case=False).mean() * 100
                else:
                    ratio = 0
                row_data.append(ratio)
            heatmap_data.append(row_data)
            
        df_heatmap = pd.DataFrame(heatmap_data, index=keywords, columns=clusters)
        
        fig_h, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(df_heatmap, cmap='YlGnBu', annot=True, fmt=".1f", ax=ax)
        ax.set_title("Keyword Frequency by Cluster (%)")
        ax.set_xlabel("Cluster ID")
        st.pyplot(fig_h)
    else:
        st.warning("数据中缺少 Cluster 列，无法绘制热图。")

# 底部
st.markdown("---")
st.caption("BioGraph App v2.1 | Powered by Streamlit & UniProt")
