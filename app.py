import streamlit as st
import pandas as pd
import networkx as nx
import os
from pyvis.network import Network
import streamlit.components.v1 as components
import gdown 

# ==========================================
# ⚙️ 页面配置与缓存
# ==========================================
st.set_page_config(page_title="BioGraph 整合分析平台", layout="wide", page_icon="🧬")

# 路径设置
# 🔴 修改点 1: 去掉 .gz 后缀，因为你 Drive 上存的可能是普通 CSV
DATA_FILE = "data/master_graph_data.csv" 
CACHE_DIR = "checkpoints"
os.makedirs(CACHE_DIR, exist_ok=True)

@st.cache_resource
def load_graph_data():
    """
    核心数据加载函数：只在启动时运行一次，之后会缓存。
    """
    # 确保 data 文件夹存在
    os.makedirs(os.path.dirname(DATA_FILE), exist_ok=True)

    # 🔴 修改点 2: 检查文件是否存在，如果文件太小(说明下载失败)也重新下载
    if not os.path.exists(DATA_FILE) or os.path.getsize(DATA_FILE) < 1000:
        # 你的 Google Drive 文件 ID
        file_id = '1tM-n8EOVCPNLg9j9JfoIgg-SFYSML5XM' 
        url = f'https://drive.google.com/uc?id={file_id}'
        
        # 显示下载提示
        with st.spinner('正在从云端下载数据 (首次运行可能需要 1-2 分钟)...'):
            # fuzzy=True 可以防止因为 Google Drive 病毒扫描提示导致的下载失败
            gdown.download(url, DATA_FILE, quiet=False, fuzzy=True)

    # 2. 读取数据
    with st.spinner('正在加载核心图谱数据...'):
        # 再次检查
        if not os.path.exists(DATA_FILE):
            st.error("下载失败。")
            return None, None

        # 🔴 修改点 3: 明确告诉 pandas 这是普通 CSV，不需要解压
        try:
            df = pd.read_csv(DATA_FILE)
        except Exception as e:
            st.error(f"读取 CSV 出错: {e}")
            # 如果读取失败，可能是文件坏了，尝试删除以便下次重启重下
            os.remove(DATA_FILE)
            return None, None
        
        # 构建 NetworkX 图
        G = nx.Graph()
        for _, row in df.iterrows():
            s, t = row['Source'], row['Target']
            # 读取属性
            score_swmed = row.get('score_Swmed', 0.0)
            is_known = row.get('in_BioGRID', False)
            loc_s = row.get('Source_Loc', 'Unknown')
            dom_s = row.get('Source_Domain', '')
            
            # 添加节点属性
            G.add_node(s, loc=loc_s, dom=dom_s)
            G.add_node(t, loc=row.get('Target_Loc', 'Unknown'), dom=row.get('Target_Domain', ''))
            
            # 添加边属性
            G.add_edge(s, t, swmed=score_swmed, known=is_known)
            
            # 计算阻力 (用于路径分析)
            resistance = 1.0 / max(score_swmed, 0.001)
            G[s][t]['resistance'] = resistance

    return G, df

# 加载数据
G, master_df = load_graph_data()

if G is None:
    st.error(f"❌ 数据加载失败。请尝试点击右上角菜单 'Clear cache' 并重启 App。")
    st.stop()

# 获取所有基因列表供下拉框使用
ALL_GENES = sorted(list(G.nodes()))

# ==========================================
# 🎨 侧边栏导航
# ==========================================
st.sidebar.title("🧬 BioGraph 分析平台")
st.sidebar.info(f"当前图谱包含:\n- 节点: {len(G.nodes)}\n- 边: {len(G.edges)}")

module = st.sidebar.radio(
    "选择功能模块:",
    [
        "🕵️‍♂️ 深度蛋白侦探",
        "🗺️ 定向路径挖掘",
        "⚖️ 全景互作分析",
        "♟️ 战略阻断模拟",
        "💧 IDR/LLPS 分析",
        "🕸️ 可视化图谱"
    ]
)

# ==========================================
# 🧩 模块 1: 深度蛋白侦探
# ==========================================
if module == "🕵️‍♂️ 深度蛋白侦探":
    st.title("🕵️‍♂️ 深度蛋白侦探 (Deep Detective)")
    st.markdown("查询目标蛋白的结合伙伴、定位与潜在机制。")

    col1, col2 = st.columns(2)
    with col1:
        default_gene = 'TP53' if 'TP53' in ALL_GENES else ALL_GENES[0] if ALL_GENES else ""
        target = st.selectbox("选择目标蛋白:", ALL_GENES, index=ALL_GENES.index(default_gene) if default_gene in ALL_GENES else 0)
    with col2:
        min_score = st.slider("最低 AI 分数阈值:", 0.0, 1.0, 0.5, 0.05)

    if st.button("开始侦查") and target:
        neighbors = []
        for n in G.neighbors(target):
            edge = G[target][n]
            score = edge.get('swmed', 0)
            is_known = edge.get('known', False)
            
            if score >= min_score or is_known:
                neighbors.append({
                    "Partner": n,
                    "Location": G.nodes[n].get('loc', 'Unk'),
                    "Type": "✅ 已知" if is_known else "🤖 预测",
                    "AI_Score": score
                })
        
        res_df = pd.DataFrame(neighbors).sort_values("AI_Score", ascending=False)
        st.write(f"### 找到 {len(res_df)} 个互作伙伴")
        st.dataframe(res_df, use_container_width=True)

# ==========================================
# 🧩 模块 2: 定向路径挖掘
# ==========================================
elif module == "🗺️ 定向路径挖掘":
    st.title("🗺️ 定向路径挖掘机")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        start_gene = st.selectbox("起点 (Start):", ALL_GENES, index=0 if ALL_GENES else 0)
    with col2:
        end_gene = st.selectbox("终点 (End):", ALL_GENES, index=1 if len(ALL_GENES)>1 else 0)
    with col3:
        via_gene = st.multiselect("必经点 (可选):", ALL_GENES)

    if st.button("挖掘路径"):
        try:
            path = []
            if not via_gene:
                # 直接路径
                path = nx.shortest_path(G, source=start_gene, target=end_gene, weight='resistance')
            else:
                # 简单处理第一个必经点
                v = via_gene[0]
                p1 = nx.shortest_path(G, start_gene, v, weight='resistance')
                p2 = nx.shortest_path(G, v, end_gene, weight='resistance')
                path = p1 + p2[1:]
            
            st.success("✅ 路径寻找成功！")
            st.code(" -> ".join(path))
            
            # 路径详细信息
            details = []
            for i in range(len(path)-1):
                u, v = path[i], path[i+1]
                edge = G[u][v]
                details.append({
                    "Step": f"{u} -> {v}",
                    "Type": "✅ 已知" if edge.get('known') else "🤖 预测",
                    "Score": edge.get('swmed', 0)
                })
            st.dataframe(pd.DataFrame(details))
            
        except nx.NetworkXNoPath:
            st.error("❌ 两点之间无通路。")
        except Exception as e:
            st.error(f"发生错误: {e}")

# ==========================================
# 🧩 模块 3: 全景互作分析
# ==========================================
elif module == "⚖️ 全景互作分析":
    st.title("⚖️ 全景互作分析 (双轨分析)")
    st.info("此处需要你在 app 中定义功能模块字典 (THEME)，为了演示，这里使用简单的定位筛选作为模拟。")
    
    locs = sorted(list(set([d.get('loc', 'Unknown') for n, d in G.nodes(data=True)])))
    
    col1, col2 = st.columns(2)
    with col1:
        loc_a = st.selectbox("区域 A (例如 Nucleus):", locs)
    with col2:
        loc_b = st.selectbox("区域 B (例如 Mitochondrion):", locs)
        
    if st.button("分析交互"):
        # 简单的模拟逻辑
        nodes_a = [n for n, d in G.nodes(data=True) if str(d.get('loc')) == loc_a]
        nodes_b = [n for n, d in G.nodes(data=True) if str(d.get('loc')) == loc_b]
        
        edges = []
        for u in nodes_a:
            for v in nodes_b:
                if G.has_edge(u, v):
                    edge = G[u][v]
                    edges.append({
                        "Source": u, "Target": v,
                        "Type": "✅ 已知" if edge.get('known') else "🤖 预测",
                        "Score": edge.get('swmed', 0)
                    })
        
        df = pd.DataFrame(edges).sort_values("Score", ascending=False)
        st.write(f"### {loc_a} <--> {loc_b} 交互列表")
        st.dataframe(df.head(50), use_container_width=True)

# ==========================================
# 🧩 模块 4: 战略阻断模拟
# ==========================================
elif module == "♟️ 战略阻断模拟":
    st.title("♟️ 战略阻断模拟 (Knockout Simulation)")
    st.markdown("模拟敲除路径上的节点，分析信号是否会被阻断或发生代偿。")
    
    col1, col2 = st.columns(2)
    with col1:
        s_node = st.selectbox("起点:", ALL_GENES, key='strat_s')
    with col2:
        e_node = st.selectbox("终点:", ALL_GENES, key='strat_e')
        
    if st.button("推演代偿机制"):
        try:
            base_path = nx.shortest_path(G, s_node, e_node, weight='resistance')
            base_cost = nx.shortest_path_length(G, s_node, e_node, weight='resistance')
            
            results = []
            targets = base_path[1:-1] # 中间节点
            
            progress_bar = st.progress(0)
            for i, target in enumerate(targets):
                G_temp = G.copy()
                G_temp.remove_node(target)
                
                try:
                    new_cost = nx.shortest_path_length(G_temp, s_node, e_node, weight='resistance')
                    impact = (new_cost - base_cost) / base_cost * 100
                    status = f"代偿代价 +{impact:.1f}%"
                except nx.NetworkXNoPath:
                    impact = 9999
                    status = "⛔ 彻底阻断"
                
                results.append({"Target": target, "Impact_Score": impact, "Status": status})
                progress_bar.progress((i+1)/len(targets))
            
            res_df = pd.DataFrame(results).sort_values("Impact_Score", ascending=False)
            
            st.write("### 靶点价值排行榜")
            st.dataframe(res_df.style.apply(lambda x: ['background-color: #ffcccc' if x['Impact_Score'] > 1000 else '' for i in x], axis=1))
            
        except Exception as e:
            st.error(f"分析出错: {e}")

# ==========================================
# 🧩 模块 5: IDR/LLPS 分析
# ==========================================
elif module == "💧 IDR/LLPS 分析":
    st.title("💧 IDR & 相分离分析仪")
    gene_input = st.text_input("输入基因名:", "FUS")
    
    if st.button("计算 IDR"):
        try:
            import metapredict as meta
            import requests
            
            with st.spinner("正在获取序列并计算..."):
                # 获取序列
                url = "https://rest.uniprot.org/uniprotkb/search"
                params = {'query': f"gene_exact:{gene_input} AND reviewed:true AND organism_id:9606", 'format': 'json', 'size': 1}
                r = requests.get(url, params=params).json()
                
                if r['results']:
                    seq = r['results'][0]['sequence']['value']
                    scores = meta.predict_disorder(seq)
                    idr_frac = sum(1 for s in scores if s > 0.5) / len(seq)
                    
                    st.metric("IDR 比例 (无序度)", f"{idr_frac:.2%}")
                    
                    if idr_frac > 0.3:
                        st.success("🌊 该蛋白具有高度无序特征，具备相分离 (LLPS) 潜力。")
                    else:
                        st.info("🪨 该蛋白结构较为刚性。")
                        
                    st.line_chart(scores)
                else:
                    st.error("未找到该基因的 UniProt 数据。")
                    
        except ImportError:
            st.warning("请在 requirements.txt 中添加 metapredict。")
        except Exception as e:
            st.error(f"计算出错: {e}")

# ==========================================
# 🧩 模块 6: 可视化图谱 (Pyvis)
# ==========================================
elif module == "🕸️ 可视化图谱":
    st.title("🕸️ 交互式子图可视化")
    center_node = st.selectbox("中心节点:", ALL_GENES)
    
    if st.button("生成图谱"):
        neighbors = list(G.neighbors(center_node))
        if len(neighbors) > 50:
            neighbors = neighbors[:50]
            st.warning("邻居过多，仅展示前 50 个。")
            
        sub_nodes = neighbors + [center_node]
        subG = G.subgraph(sub_nodes)
        
        net = Network(height="600px", width="100%", bgcolor="#222222", font_color="white")
        
        for n in subG.nodes():
            color = "#ff4757" if n == center_node else "#1e90ff"
            net.add_node(n, label=n, color=color, title=f"Loc: {subG.nodes[n].get('loc')}")
            
        for u, v in subG.edges():
            edge = subG[u][v]
            color = "#7bed9f" if edge.get('known') else "#ffa502"
            net.add_edge(u, v, color=color, width=2 if edge.get('known') else 1)
            
        tmp_path = "tmp_network.html"
        net.save_graph(tmp_path)
        
        with open(tmp_path, 'r', encoding='utf-8') as f:
            source_code = f.read()
            components.html(source_code, height=620)
