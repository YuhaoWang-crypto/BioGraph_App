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
# 1. 页面配置
# ==========================================
st.set_page_config(
    page_title="BioGraph v3.0: Protein Lifecycle Explorer",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🧬 BioGraph v3.0: 多维蛋白质组学流形分析平台")
st.markdown("""
**真实数据驱动**：集成 UMAP 流形投影、功能网络拓扑与深度特征挖掘。
""")

# ==========================================
# 2. 数据加载与预处理 (核心升级)
# ==========================================
@st.cache_data
def load_data():
    try:
        # 读取压缩文件
        df = pd.read_csv("final_analysis_result.csv.gz", compression='gzip')
        
        # --- 基础清洗 ---
        if 'cc_function' in df.columns:
            df['cc_function'] = df['cc_function'].fillna('Unknown')
        if 'Gene_Symbol' in df.columns:
            df['Gene_Symbol'] = df['Gene_Symbol'].fillna('Unknown')
        if 'Real_Protein_HalfLife_Hours' in df.columns:
            df['Real_Protein_HalfLife_Hours'] = df['Real_Protein_HalfLife_Hours'].fillna(0)
            
        # --- 生成辅助标签 1: 稳定性 ---
        if 'Real_Protein_HalfLife_Hours' in df.columns:
            df['Stability_Level'] = pd.cut(df['Real_Protein_HalfLife_Hours'], 
                                           bins=[-1, 10, 50, 10000], 
                                           labels=['Short (<10h)', 'Medium', 'Long (>50h)']).astype(str)
                                           
        # --- 生成辅助标签 2: 癌症相关性 ---
        if 'Is_Cancer' not in df.columns and 'cc_function' in df.columns:
             df['Is_Cancer'] = df['cc_function'].str.contains('cancer|tumor', case=False).map({True:'Yes', False:'No'})

        # --- 生成辅助标签 3: 细胞位置 (Auto_Location) ---
        # 即使 CSV 里没有，这里也会自动算出来，保证 Tab 1 能用
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
        
        return df
    except FileNotFoundError:
        st.error("❌ 未找到数据文件！请上传 final_analysis_result.csv.gz")
        return pd.DataFrame()
    except Exception as e:
        st.error(f"❌ 数据读取错误: {e}")
        return pd.DataFrame()

df_main = load_data()
if df_main.empty: st.stop()

# ==========================================
# 3. 计算 PCA Loadings
# ==========================================
@st.cache_data
def calculate_pca_loadings(df):
    features = ['Real_Protein_HalfLife_Hours', 'mRNA_Expression', 'Circadian_Amplitude']
    valid_features = [f for f in features if f in df.columns]
    if len(valid_features) < 2: return None, None
    X = np.log1p(df[valid_features].fillna(0))
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    pca = PCA(n_components=2)
    pca.fit(X_scaled)
    loadings = pd.DataFrame(pca.components_.T, columns=['PCA1', 'PCA2'], index=valid_features)
    return loadings, valid_features

df_loadings, used_features = calculate_pca_loadings(df_main)

# ==========================================
# 4. 侧边栏设置
# ==========================================
st.sidebar.header("🔍 全局控制")

# Cluster 筛选
if 'Cluster' in df_main.columns:
    all_clusters = sorted(df_main['Cluster'].unique())
    selected_clusters = st.sidebar.multiselect("筛选 Cluster", all_clusters, default=all_clusters)
    df_filtered = df_main[df_main['Cluster'].isin(selected_clusters)]
else:
    df_filtered = df_main

# 侧边栏搜索 (联动 Tab 1 高亮)
sidebar_search = st.sidebar.text_input("在全景图中高亮基因 (如 TP53)", "").upper()
st.sidebar.info(f"展示: {len(df_filtered)} 条")

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

# --- Tab 1: 全景流形 (已添加 Auto_Location) ---
with tab1:
    col1, col2 = st.columns([1, 4])
    with col1:
        # 动态检测可用列，确保 Auto_Location 在里面
        options = ['Cluster', 'Stability_Level', 'Auto_Location', 'Is_Cancer', 'N_Term_AA']
        valid_options = [c for c in options if c in df_filtered.columns]
        
        color_col = st.radio("上色依据 (Color By):", valid_options, index=0)

    with col2:
        if 'UMAP_X' in df_filtered.columns:
            fig = px.scatter(
                df_filtered, 
                x='UMAP_X', y='UMAP_Y', 
                color=color_col,
                hover_data=['Gene_Symbol', 'Real_Protein_HalfLife_Hours', 'Auto_Location'],
                title=f"Functional Manifold (Colored by {color_col})",
                height=650,
                template="plotly_white",
                opacity=0.6,
                color_discrete_sequence=px.colors.qualitative.Bold
            )
            
            # 侧边栏搜索高亮
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
            st.warning("缺少 UMAP 坐标数据。")

# --- Tab 2: 详情雷达 (已添加基因选择框) ---
with tab2:
    st.markdown("### 单蛋白详细档案查询")
    
    # 1. 基因选择器 (下拉框 + 搜索)
    all_genes = sorted(df_main['Gene_Symbol'].unique())
    # 智能联动：如果侧边栏搜了，这里默认选中它
    default_idx = 0
    if sidebar_search and sidebar_search in all_genes:
        default_idx = all_genes.index(sidebar_search)
        
    selected_gene = st.selectbox("选择或输入基因名:", all_genes, index=default_idx)
    
    # 2. 展示详情
    if selected_gene:
        row = df_main[df_main['Gene_Symbol'] == selected_gene].iloc[0]
        
        c1, c2 = st.columns([1, 1])
        with c1:
            # 防弹版静态定位图
            fig_loc, ax = plt.subplots(figsize=(8, 6))
            # 背景点
            sns.scatterplot(data=df_main.sample(min(2000, len(df_main))), x='PCA1', y='PCA2', 
                            color='lightgrey', s=10, alpha=0.3, ax=ax, label='Background')
            # 选中点
            ax.scatter(row['PCA1'], row['PCA2'], color='red', s=200, marker='*', edgecolors='black', zorder=10)
            ax.text(row['PCA1'], row['PCA2']+0.3, selected_gene, color='red', fontweight='bold', ha='center')
            ax.set_title(f"Protein Position in PCA Space")
            ax.set_xlabel("PCA1 (Stability)")
            ax.set_ylabel("PCA2 (Dynamics)")
            st.pyplot(fig_loc)
            
        with c2:
            st.subheader(f"🧬 {selected_gene}")
            st.success(f"Cluster: {row.get('Cluster', 'N/A')}")
            st.info(f"Location: {row.get('Auto_Location', 'Unknown')}")
            
            st.metric("真实半衰期", f"{row.get('Real_Protein_HalfLife_Hours', 0):.1f} h")
            st.metric("N端氨基酸", f"{row.get('N_Term_AA', 'N/A')}")
            
            st.markdown("**功能描述:**")
            st.caption(row.get('cc_function', 'No description.'))

# --- Tab 3: 功能网络 (已改为下拉菜单) ---
with tab3:
    st.markdown("### 功能模块共现网络 (Co-occurrence Network)")
    
    # 预定义模块列表 (用户无需手打)
    modules = [
        'Mitochondria (线粒体)', 
        'Nucleus (细胞核)', 
        'Plasma Membrane (细胞膜)', 
        'Ribosome (核糖体)',
        'Cytoskeleton (细胞骨架)',
        'Kinase (激酶)',
        'Ubiquitin (泛素)',
        'DNA Repair (DNA修复)',
        'Cell Cycle (细胞周期)',
        'Apoptosis (凋亡)',
        'Immune Response (免疫)'
    ]
    
    selected_module_label = st.selectbox("选择感兴趣的功能模块:", modules)
    
    # 提取关键词 (去除括号里的中文)
    keyword = selected_module_label.split(' (')[0]
    
    if keyword:
        # 筛选数据
        subset = df_main[df_main['cc_function'].str.contains(keyword, case=False, na=False)].head(100)
        
        if len(subset) > 2:
            G = nx.Graph()
            genes = subset['Gene_Symbol'].tolist()
            hls = subset['Real_Protein_HalfLife_Hours'].tolist()
            
            # 建图逻辑：半衰期接近的连线
            for i in range(len(genes)):
                G.add_node(genes[i])
                for j in range(i+1, len(genes)):
                    if abs(hls[i] - hls[j]) < 2.0: # 连线阈值
                        G.add_edge(genes[i], genes[j])
            
            # 绘图
            fig_net, ax = plt.subplots(figsize=(10, 7))
            pos = nx.spring_layout(G, k=0.2, seed=42)
            
            # 节点颜色映射半衰期
            nodes = nx.draw_networkx_nodes(G, pos, node_size=50, 
                                         node_color=hls, cmap='viridis', 
                                         alpha=0.8, ax=ax)
            nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color='gray', ax=ax)
            nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
            
            plt.colorbar(nodes, label='Half-life (Hours)')
            ax.set_title(f"Network: {keyword} (Color = Half-life)")
            ax.axis('off')
            st.pyplot(fig_net)
            
            st.markdown(f"**节点数:** {len(genes)} | **连线逻辑:** 半衰期差异 < 2小时")
        else:
            st.warning(f"模块 '{keyword}' 中的蛋白数量过少 (<3)，无法构建网络。")

# --- Tab 4 & 5: 保持原样 (PCA Loadings & Heatmap) ---
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
            st.pyplot(fig_l)

with tab5:
    st.write("### 关键词富集热图")
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
                if len(sub)>0: ratio = sub['cc_function'].str.contains(k, case=False).mean()*100
                else: ratio = 0
            row_data.append(row_data)
        
        df_heatmap = pd.DataFrame(heatmap_data, index=keywords, columns=clusters)
        fig_h, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(df_heatmap, cmap='YlGnBu', annot=True, fmt=".1f", ax=ax)
        st.pyplot(fig_h)

st.markdown("---")
st.caption("BioGraph v3.0 | Streamlit")
