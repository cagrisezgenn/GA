# 1) MATLAB içeren resmi imaj (kampüs lisansı için uygun)
FROM mathworks/matlab:r2024a

# 2) Faydalı araçlar
RUN apt-get update && apt-get install -y --no-install-recommends \
    git git-lfs curl ca-certificates tmux && \
    rm -rf /var/lib/apt/lists/*
