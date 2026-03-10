FROM continuumio/miniconda3:latest

# Install ViennaRNA (for RNAfold) and Perl (for miRCheck)
RUN conda install -c conda-forge -y viennarna=2.5.1 perl && \
    conda clean -afy

# Install Python dependencies
RUN pip install --no-cache-dir \
    numpy \
    pandas \
    scipy \
    statsmodels \
    scikit-learn \
    matplotlib

# Copy siWalk
COPY . /siWalk
WORKDIR /siWalk/src

# Make bundled binaries executable
RUN chmod +x /siWalk/lib/miranda

ENTRYPOINT ["python", "siWalk_predict_siRNA_location.py"]
