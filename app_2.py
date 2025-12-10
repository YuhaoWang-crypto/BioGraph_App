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

# ==========================================
# 1. 页面基础配置
# ==========================================
st.set_page_config(
    page_title="BioGraph v4.0: Stable Release",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🧬 BioGraph v4.0: 蛋白质组学全景分析平台")
st.markdown("""
**稳定版更新说明**：
1. 修复了网络图谱 (Tab 3) 中因基因重复导致的绘图报错。
2. 增强了数据加载的容错性，自动补全缺失的分类标签。
""")

# ==========================================
# 2. 数据加载与智能预处理 (核心引擎)
# ==========================================
@st.cache_data
def load_data():
    try:
        # 读取压缩数据
        df = pd.read_csv("final_analysis_result.csv.gz", compression='gzip')
        
        # --- A. 基础清洗 (防止空值报错) ---
        str_cols = ['cc_function', 'Gene_Symbol', 'N_Term_AA', 'Processing_Type']
        for col in str_cols:
            if col in df.columns:
                df[col] = df[col].fillna('Unknown')
        
        if 'Real_Protein_HalfLife_Hours' in df.columns:
            df['Real_Protein_HalfLife_Hours'] = df['Real_Protein_HalfLife_Hours'].fillna(0)
            
        # --- B. 智能标签补全 (如果CSV里没有这些列，现场算) ---
        
        # 1. 细胞位置 (Auto_Location)
        if 'Auto_Location' not in df.columns and 'cc_function' in df.columns:
            def get_loc(text):
                t = str(text).lower()
                if 'mitoch' in t: return 'Mitochondria'
                if 'nucleus' in t or 'nuclear' in t: return 'Nucleus'
                if 'membrane' in t and 'plasma' in t: return 'Plasma Membrane'
                if 'ribosom' in t: return 'Ribosome'
                if 'endoplasmic' in t or 'reticulum' in t: return 'ER'
                if 'golgi' in t: return 'Golgi'
                if 'secreted' in t: return 'Secreted'
                return 'Cytoplasm/Other'
            df['Auto_Location'] = df['cc_function'].apply(get_loc)

        # 2. 癌症相关性 (Is_Cancer)
        if 'Is_Cancer' not in df.columns and 'cc_function' in df.columns:
             df['Is_Cancer'] = df['cc_function'].str.contains('cancer|tumor', case=False).map({True:'Yes', False:'No'})

        # 3. 稳定性分级 (Stability_Level)
        if 'Stability_Level' not in df.columns and 'Real_Protein_HalfLife_Hours' in df.columns:
            df['Stability_Level'] = pd.cut(df['Real_Protein_HalfLife_Hours'], 
                                           bins=[-1, 10, 50, 100000], 
                                           labels=['Short (<10h)', 'Medium', 'Long (>50h)']).astype(str)

        return df
        
    except FileNotFoundError:
        st.error("❌ 严重错误：未找到数据文件 `final_analysis_result.csv.gz`。请确保文件已上传到 GitHub 仓库根目录。")
        return pd.DataFrame()
    except Exception as e:
        st.error(f"❌ 数据读取未知错误: {e}")
        return pd.DataFrame()

# 加载数据
df_main = load_data()

# 如果数据加载失败，停止运行
if df_main.empty:
    st.stop()

# ==========================================
# 3. 辅助计算：PCA Loadings
# ==========================================
@st.cache_data
def calculate_pca_loadings(df):
    features = ['Real_Protein_HalfLife_Hours', 'mRNA_Expression', 'Circadian_Amplitude']
    valid_features = [f for f in features if f in df.columns]
    
    if len(valid_features) < 2: return None, None
    
    # 简单的 PCA 重算用于画向量图
    X = np.log1p(df[valid_features].fillna(0))
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=2)
    pca.fit(X_scaled)
    
    loadings = pd.DataFrame(pca.components_.T, columns=['PCA1', 'PCA2'], index=valid_features)
    return loadings, valid_features

df_loadings, used_features = calculate_pca_loadings(df_main)

# ==========================================
# 4. 侧边栏：全局过滤器
# ==========================================
st.sidebar.header("🔍 全局筛选")

# Cluster 筛选器
if 'Cluster' in df_main.columns:
    all_clusters = sorted(df_main['Cluster'].unique())
    selected_clusters = st.sidebar.multiselect("筛选 Cluster", all_clusters, default=all_clusters)
    df_filtered = df_main[df_main['Cluster'].isin(selected_clusters)]
else:
    df_filtered = df_main

# 全局搜索框
sidebar_search = st.sidebar.text_input("全景图基因高亮 (如 TP53)", "").upper()

st.sidebar.info(f"当前展示: {len(df_filtered)} 条数据")

# ==========================================
# 5. 主界面 Tabs
# ==========================================
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🌌 全景流形 (UMAP)", 
    "🔍 蛋白雷达 (Detail)", 
    "🕸️ 功能网络 (Network)", 
    "📉 PCA 解密",
    "🧪 专业富集"
])

# --- Tab 1: 全景流形图 ---
with tab1:
    col1, col2 = st.columns([1, 4])
    with col1:
        # 动态检测可用列
        options = ['Cluster', 'Auto_Location', 'Stability_Level', 'Is_Cancer', 'N_Term_AA']
        valid_opts = [o for o in options if o in df_filtered.columns]
        color_col = st.radio("上色依据:", valid_opts, index=0 if valid_opts else None)

    with col2:
        if 'UMAP_X' in df_filtered.columns and color_col:
            fig = px.scatter(
                df_filtered, 
                x='UMAP_X', y='UMAP_Y', 
                color=color_col,
                hover_data=['Gene_Symbol', 'Real_Protein_HalfLife_Hours', 'Auto_Location'],
                title=f"Functional Manifold (Colored by {color_col})",
                height=600,
                template="plotly_white",
                opacity=0.6,
                color_discrete_sequence=px.colors.qualitative.Bold
            )
            
            # 高亮搜索点
            if sidebar_search and not df_filtered[df_filtered['Gene_Symbol'] == sidebar_search].empty:
                row = df_filtered[df_filtered['Gene_Symbol'] == sidebar_search].iloc[0]
                fig.add_trace(go.Scatter(
                    x=[row['UMAP_X']], y=[row['UMAP_Y']],
                    mode='markers+text',
                    marker=dict(size=25, color='red', symbol='star', line=dict(width=2, color='black')),
                    text=[sidebar_search],
                    textposition="top center",
                    name='Searched'
                ))
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("缺少 UMAP 坐标数据，无法绘图。")

# --- Tab 2: 详情雷达 ---
with tab2:
    # 基因选择器逻辑：优先使用侧边栏搜索结果
    all_genes = sorted(df_main['Gene_Symbol'].unique())
    default_idx = 0
    if sidebar_search and sidebar_search in all_genes:
        default_idx = all_genes.index(sidebar_search)
        
    selected_gene = st.selectbox("选择或输入基因名查看详情:", all_genes, index=default_idx)
    
    if selected_gene:
        row = df_main[df_main['Gene_Symbol'] == selected_gene].iloc[0]
        
        c1, c2 = st.columns([1, 1])
        with c1:
            # 静态定位图 (Matplotlib) - 这种图最稳定，不会崩
            if 'PCA1' in df_main.columns:
                fig_loc, ax = plt.subplots(figsize=(8, 6))
                # 背景
                sns.scatterplot(data=df_main.sample(min(3000, len(df_main))), x='PCA1', y='PCA2', 
                                color='lightgrey', s=10, alpha=0.3, ax=ax)
                # 高亮
                ax.scatter(row['PCA1'], row['PCA2'], color='red', s=200, marker='*', edgecolors='black', zorder=10)
                ax.text(row['PCA1'], row['PCA2']+0.3, selected_gene, color='red', fontweight='bold', ha='center')
                ax.set_title("Position in PCA Space")
                ax.set_xlabel("PCA1 (Stability)")
                ax.set_ylabel("PCA2 (Dynamics)")
                st.pyplot(fig_loc)
            else:
                st.warning("缺少 PCA 数据。")
            
        with c2:
            st.subheader(f"🧬 {selected_gene}")
            st.write(f"**Cluster ID:** {row.get('Cluster', 'N/A')}")
            st.write(f"**Location:** {row.get('Auto_Location', 'N/A')}")
            st.metric("真实半衰期", f"{row.get('Real_Protein_HalfLife_Hours', 0):.1f} h")
            st.metric("mRNA 表达量", f"{row.get('mRNA_Expression', 0):.2f}")
            st.markdown("**功能描述:**")
            st.info(row.get('cc_function', 'No description.'))

# --- Tab 3: 功能网络 (关键修复版) ---
with tab3:
    st.markdown("### 功能模块共现网络 (Co-occurrence Network)")
    
    modules = [
        'Mitochondria (线粒体)', 'Nucleus (细胞核)', 'Plasma Membrane (细胞膜)', 
        'Ribosome (核糖体)', 'Cytoskeleton (细胞骨架)', 'Kinase (激酶)',
        'Ubiquitin (泛素)', 'DNA Repair (DNA修复)', 'Cell Cycle (细胞周期)',
        'Apoptosis (凋亡)', 'Immune Response (免疫)'
    ]
    
    selected_module = st.selectbox("选择功能模块:", modules)
    keyword = selected_module.split(' (')[0] # 提取英文关键词
    
    if keyword:
        # 1. 筛选数据
        subset = df_main[df_main['cc_function'].str.contains(keyword, case=False, na=False)]
        
        # === 核心修复：去重 ===
        # 确保每个基因只出现一次，防止绘图时节点与颜色数量不匹配
        subset = subset.drop_duplicates(subset=['Gene_Symbol'])
        
        # 限制数量防止浏览器卡顿
        subset = subset.head(80)
        
        if len(subset) > 2:
            G = nx.Graph()
            genes = subset['Gene_Symbol'].tolist()
            hls = subset['Real_Protein_HalfLife_Hours'].tolist()
            
            # 2. 建图 (将属性写入节点)
            for i in range(len(genes)):
                G.add_node(genes[i], half_life=hls[i]) # 存属性！
                for j in range(i+1, len(genes)):
                    if abs(hls[i] - hls[j]) < 2.0:
                        G.add_edge(genes[i], genes[j])
            
            # 3. 准备绘图 (从 Graph 对象中按顺序提取颜色)
            # 这样能保证颜色列表与 G.nodes() 的顺序绝对一致
            node_colors = [G.nodes[n]['half_life'] for n in G.nodes()]
            
            # 4. 绘图
            fig_net, ax = plt.subplots(figsize=(10, 7))
            pos = nx.spring_layout(G, k=0.2, seed=42)
            
            nodes = nx.draw_networkx_nodes(G, pos, 
                                         node_size=60, 
                                         node_color=node_colors, 
                                         cmap='viridis', 
                                         alpha=0.8, ax=ax)
            
            nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color='gray', ax=ax)
            
            # 节点少的时候才显示标签
            if len(G.nodes()) < 60:
                nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
            
            # 添加 Colorbar
            if len(node_colors) > 0:
                plt.colorbar(nodes, label='Half-life (Hours)', ax=ax)
                
            ax.set_title(f"Network: {keyword} (Color = Half-life)")
            ax.axis('off')
            st.pyplot(fig_net)
            
            st.caption(f"展示节点数: {len(genes)} | 连线规则: 半衰期差异 < 2h")
        else:
            st.warning(f"模块 '{keyword}' 数据过少 (<3)，无法构建网络。")

# --- Tab 4: PCA Loadings ---
with tab4:
    if df_loadings is not None:
        col_l, col_r = st.columns(2)
        with col_l:
            st.dataframe(df_loadings.style.background_gradient(cmap='RdBu'))
        with col_r:
            fig_l, ax = plt.subplots(figsize=(6, 6))
            ax.axhline(0, color='grey', linestyle='--')
            ax.axvline(0, color='grey', linestyle='--')
            for i, feat in enumerate(df_loadings.index):
                ax.arrow(0, 0, df_loadings.iloc[i, 0], df_loadings.iloc[i, 1], color='red', width=0.01)
                ax.text(df_loadings.iloc[i, 0]*1.2, df_loadings.iloc[i, 1]*1.2, feat)
            ax.set_xlabel("PCA 1")
            ax.set_ylabel("PCA 2")
            st.pyplot(fig_l)
    else:
        st.warning("无法计算 PCA Loadings，可能缺少数值列。")

# --- Tab 5: 富集热图 ---
with tab5:
    st.write("### 关键词富集热图 (Keyword Enrichment)")
    
    BIO_DICT = {
        'Loc': ['mitochondrion', 'nucleus', 'membrane', 'secreted'],
        'Func': ['kinase', 'transcription', 'transport', 'metabolism'],
        'Struct': ['zinc', 'finger', 'domain']
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
        ax.set_xlabel("Cluster ID")
        st.pyplot(fig_h)
    else:
        st.warning("缺少 Cluster 列，无法绘制热图。")

# 底部版权
st.markdown("---")
st.caption("BioGraph v4.0 Stable | Powered by Streamlit")
