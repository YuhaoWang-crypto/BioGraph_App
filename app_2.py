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
    page_title="BioGraph v4.5: Stable Interactive",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🧬 BioGraph v4.5: 交互式蛋白质组学全景平台")
st.markdown("""
**稳定增强版**：基于 v4.0 内核。
修复了网络图报错，并集成了交互式网络 (Plotly) 与动态热图分析。
""")

# ==========================================
# 2. 数据加载与预处理 (沿用 v4.0 的高容错逻辑)
# ==========================================
@st.cache_data
def load_data():
    try:
        df = pd.read_csv("final_analysis_result.csv.gz", compression='gzip')
        
        # 基础清洗
        str_cols = ['cc_function', 'Gene_Symbol', 'N_Term_AA', 'Processing_Type']
        for col in str_cols:
            if col in df.columns:
                df[col] = df[col].fillna('Unknown')
        
        if 'Real_Protein_HalfLife_Hours' in df.columns:
            df['Real_Protein_HalfLife_Hours'] = df['Real_Protein_HalfLife_Hours'].fillna(0)
            
        # 标签补全
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

        if 'Is_Cancer' not in df.columns and 'cc_function' in df.columns:
             df['Is_Cancer'] = df['cc_function'].str.contains('cancer|tumor', case=False).map({True:'Yes', False:'No'})

        if 'Stability_Level' not in df.columns and 'Real_Protein_HalfLife_Hours' in df.columns:
            df['Stability_Level'] = pd.cut(df['Real_Protein_HalfLife_Hours'], 
                                           bins=[-1, 10, 50, 100000], 
                                           labels=['Short (<10h)', 'Medium', 'Long (>50h)']).astype(str)

        return df
    except Exception as e:
        st.error(f"❌ 数据读取错误: {e}")
        return pd.DataFrame()

df_main = load_data()
if df_main.empty: st.stop()

# ==========================================
# 3. PCA Loadings 计算
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
# 4. 侧边栏
# ==========================================
st.sidebar.header("🔍 全局筛选")
if 'Cluster' in df_main.columns:
    all_clusters = sorted(df_main['Cluster'].unique())
    selected_clusters = st.sidebar.multiselect("筛选 Cluster", all_clusters, default=all_clusters)
    df_filtered = df_main[df_main['Cluster'].isin(selected_clusters)]
else:
    df_filtered = df_main

sidebar_search = st.sidebar.text_input("全景图基因高亮 (如 TP53)", "").upper()
st.sidebar.info(f"展示: {len(df_filtered)} 条数据")

# ==========================================
# 5. 主界面 Tabs
# ==========================================
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🌌 全景流形 (UMAP)", 
    "🔍 蛋白雷达 (Detail)", 
    "🕸️ 交互网络 (Interactive)", 
    "📉 PCA 解密",
    "🧪 动态富集 (Heatmap)"
])

# --- Tab 1: 全景流形 (保持 v4 的 Plotly 实现) ---
with tab1:
    col1, col2 = st.columns([1, 4])
    with col1:
        options = ['Cluster', 'Auto_Location', 'Stability_Level', 'Is_Cancer', 'N_Term_AA']
        valid_opts = [o for o in options if o in df_filtered.columns]
        color_col = st.radio("上色依据:", valid_opts, index=0 if valid_opts else None)

    with col2:
        if 'UMAP_X' in df_filtered.columns and color_col:
            fig = px.scatter(
                df_filtered, x='UMAP_X', y='UMAP_Y', color=color_col,
                hover_data=['Gene_Symbol', 'Real_Protein_HalfLife_Hours', 'Auto_Location'],
                title=f"Functional Manifold (Colored by {color_col})",
                height=600, template="plotly_white", opacity=0.6,
                color_discrete_sequence=px.colors.qualitative.Bold
            )
            if sidebar_search and not df_filtered[df_filtered['Gene_Symbol'] == sidebar_search].empty:
                row = df_filtered[df_filtered['Gene_Symbol'] == sidebar_search].iloc[0]
                fig.add_trace(go.Scatter(
                    x=[row['UMAP_X']], y=[row['UMAP_Y']], mode='markers+text',
                    marker=dict(size=25, color='red', symbol='star', line=dict(width=2, color='black')),
                    text=[sidebar_search], textposition="top center", name='Searched'
                ))
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("缺少 UMAP 坐标。")

# --- Tab 2: 详情雷达 (保持 v4 的静态图，最稳定) ---
with tab2:
    all_genes = sorted(df_main['Gene_Symbol'].unique())
    default_idx = all_genes.index(sidebar_search) if sidebar_search and sidebar_search in all_genes else 0
    selected_gene = st.selectbox("选择或输入基因名查看详情:", all_genes, index=default_idx)
    
    if selected_gene:
        row = df_main[df_main['Gene_Symbol'] == selected_gene].iloc[0]
        c1, c2 = st.columns([1, 1])
        with c1:
            if 'PCA1' in df_main.columns:
                fig_loc, ax = plt.subplots(figsize=(8, 6))
                sns.scatterplot(data=df_main.sample(min(3000, len(df_main))), x='PCA1', y='PCA2', 
                                color='lightgrey', s=10, alpha=0.3, ax=ax)
                ax.scatter(row['PCA1'], row['PCA2'], color='red', s=200, marker='*', edgecolors='black', zorder=10)
                ax.text(row['PCA1'], row['PCA2']+0.3, selected_gene, color='red', fontweight='bold', ha='center')
                ax.set_title("Position in PCA Space")
                st.pyplot(fig_loc)
        with c2:
            st.subheader(f"🧬 {selected_gene}")
            st.write(f"**Cluster:** {row.get('Cluster', 'N/A')} | **Loc:** {row.get('Auto_Location', 'N/A')}")
            st.metric("真实半衰期", f"{row.get('Real_Protein_HalfLife_Hours', 0):.1f} h")
            st.info(row.get('cc_function', 'No description.'))

# --- Tab 3: 交互式网络 (修复 KeyError 版) ---
with tab3:
    st.markdown("### 🕸️ 交互式功能共现网络")
    st.caption("基于半衰期相似性构建网络。支持缩放、悬停查看详情。")
    
    modules = [
        'Mitochondria (线粒体)', 'Nucleus (细胞核)', 'Plasma Membrane (细胞膜)', 
        'Ribosome (核糖体)', 'Cytoskeleton (细胞骨架)', 'Kinase (激酶)',
        'Ubiquitin (泛素)', 'DNA Repair (DNA修复)', 'Cell Cycle (细胞周期)',
        'Apoptosis (凋亡)', 'Immune Response (免疫)'
    ]
    selected_module = st.selectbox("选择功能模块:", modules)
    keyword = selected_module.split(' (')[0]
    
    if keyword:
        subset = df_main[df_main['cc_function'].str.contains(keyword, case=False, na=False)]
        subset = subset.drop_duplicates(subset=['Gene_Symbol']).head(80)
        
        if len(subset) > 2:
            G = nx.Graph()
            
            # === 关键修复：全部转换为 list，避免索引对不上导致的 KeyError ===
            genes_list = subset['Gene_Symbol'].tolist()
            hls_list = subset['Real_Protein_HalfLife_Hours'].tolist()
            funcs_list = subset['cc_function'].astype(str).tolist()
            
            for i in range(len(genes_list)):
                # 将属性写入节点
                G.add_node(genes_list[i], hl=hls_list[i], desc=funcs_list[i])
                for j in range(i+1, len(genes_list)):
                    if abs(hls_list[i] - hls_list[j]) < 2.0:
                        G.add_edge(genes_list[i], genes_list[j])
            
            # 布局
            pos = nx.spring_layout(G, k=0.3, seed=42)
            
            # Plotly 绘图：边
            edge_x = []
            edge_y = []
            for edge in G.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])

            edge_trace = go.Scatter(
                x=edge_x, y=edge_y,
                line=dict(width=0.5, color='#888'),
                hoverinfo='none',
                mode='lines')

            # Plotly 绘图：点
            node_x = []
            node_y = []
            node_text = []
            node_color = []
            
            for node in G.nodes():
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)
                
                # 获取属性
                hl = G.nodes[node]['hl']
                desc = G.nodes[node]['desc'][:100] + "..." # 截断长文本
                
                node_text.append(f"<b>{node}</b><br>HL: {hl:.1f}h<br>{desc}")
                node_color.append(hl)

            node_trace = go.Scatter(
                x=node_x, y=node_y,
                mode='markers',
                hoverinfo='text',
                text=node_text,
                marker=dict(
                    showscale=True,
                    colorscale='Viridis',
                    color=node_color,
                    size=15,
                    colorbar=dict(title='Half-Life (h)'),
                    line_width=1))

            fig_net = go.Figure(data=[edge_trace, node_trace],
                         layout=go.Layout(
                            title=f"Network: {keyword}",
                            showlegend=False,
                            hovermode='closest',
                            margin=dict(b=20,l=5,r=5,t=40),
                            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            height=600,
                            template='plotly_white'
                         ))
            st.plotly_chart(fig_net, use_container_width=True)
        else:
            st.warning("相关蛋白数量过少，无法构建网络。")

# --- Tab 4: PCA Loadings ---
with tab4:
    if df_loadings is not None:
        c1, c2 = st.columns(2)
        with c1: st.dataframe(df_loadings.style.background_gradient(cmap='RdBu'))
        with c2: 
            fig_l, ax = plt.subplots(figsize=(6,6))
            ax.axhline(0, c='grey', ls='--'); ax.axvline(0, c='grey', ls='--')
            for i, f in enumerate(df_loadings.index):
                ax.arrow(0,0, df_loadings.iloc[i,0], df_loadings.iloc[i,1], color='r', width=0.01)
                ax.text(df_loadings.iloc[i,0]*1.2, df_loadings.iloc[i,1]*1.2, f)
            st.pyplot(fig_l)

# --- Tab 5: 动态富集 (升级为 Plotly 交互热图) ---
with tab5:
    st.markdown("### 🧪 动态关键词富集分析")
    
    FULL_DICT = {
        'Location': ['mitochondrion', 'nucleus', 'membrane', 'cytoplasm', 'secreted', 'golgi', 'ER'],
        'Function': ['kinase', 'transcription', 'transport', 'metabolism', 'receptor', 'chaperone'],
        'Process': ['cell cycle', 'apoptosis', 'immune', 'signaling', 'dna repair']
    }
    
    selected_cats = st.multiselect("选择分析维度:", list(FULL_DICT.keys()), default=['Location', 'Function'])
    
    target_kws = []
    for cat in selected_cats:
        target_kws.extend(FULL_DICT[cat])
    
    if target_kws and 'Cluster' in df_main.columns:
        clusters = sorted(df_main['Cluster'].unique())
        heatmap_data = []
        
        for k in target_kws:
            row_data = []
            for c in clusters:
                sub = df_main[df_main['Cluster'] == c]
                ratio = sub['cc_function'].str.contains(k, case=False).mean() * 100 if len(sub)>0 else 0
                row_data.append(ratio)
            heatmap_data.append(row_data)
            
        df_hm = pd.DataFrame(heatmap_data, index=target_kws, columns=[f"Cluster {c}" for c in clusters])
        
        fig_h = px.imshow(
            df_hm,
            labels=dict(x="Cluster", y="Keyword", color="%"),
            color_continuous_scale='YlGnBu',
            aspect="auto",
            title="Keyword Enrichment (%)"
        )
        fig_h.update_traces(text=df_hm.round(1).values, texttemplate="%{text}%")
        st.plotly_chart(fig_h, use_container_width=True)

st.markdown("---")
st.caption("BioGraph v4.5 Stable | Powered by Streamlit")
