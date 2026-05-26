import numpy as np
import pickle
from pomegranate import *
from time import time
from figsmod.hmm import initialize_hmm
from figsmod.protein import Protein


def read_biogrid(biogrid_path, center):
    with open(biogrid_path, 'rb') as f:
        biogrid_net = pickle.load(f)
    input_proteins_id = []
    if len(biogrid_net[center]) <= 50 and len(biogrid_net[center]) >= 10:
        for u, v, data in biogrid_net.edges(center, data=True):
            if v == '-' or u == '-':
                continue
            input_proteins_id.append(v)
    else:
        for u, v, data in biogrid_net.edges(center, data=True):
            if v == '-' or u == '-':
                continue
            if len(data['BiogridIDs']) >= 2 or data['low_throughput']:
                input_proteins_id.append(v)
    return input_proteins_id


def load_protein_database(protein_database_path, input_proteins_id):
    with open(protein_database_path, 'rb') as f:
        protein_database = pickle.load(f)
    proteins = []
    missing = []
    for id in input_proteins_id:
        if id in protein_database:
            proteins.append(protein_database[id])
        else:
            missing.append(id)
    if missing:
        print(f"Warning: {len(missing)} protein(s) not found in database and skipped:")
        for m in missing:
            print(f"  {m}")
    if not proteins:
        raise ValueError("No valid proteins found. Check your input file and protein database.")
    return proteins


def eval_score(model, k):
    score = 0
    for state in model.states:
        if state.name.startswith('m'):
            info = np.log2(20)
            for r, p in state.distribution.parameters[0].items():
                info += p * np.log2(p, where=p > 0)
            score += info
    return score / k


def search(proteins, seed, iterations, k, T, w_plm, w_disorder, w_rlc, w_rsa):
    np.random.seed(seed)
    # compute biased sampling probability 
    for protein in proteins:
        plm = protein.plm
        disorder = protein.disorder
        rlc = protein.rlc
        rsa = protein.rsa
        weight = w_plm * plm + w_disorder * disorder + w_rlc * rlc + w_rsa * rsa
        prob = np.zeros(len(weight) - k + 1)
        for i in range(len(weight) - k + 1):
            prob[i] = np.sum(weight[i:i+k])
        prob = prob - min(prob) + 1e-10   # shifting the distribution to avoid non-positive
        prob = prob / np.sum(prob)
        protein.prob = prob


    # initialize positions
    positions = [0] * len(proteins)
    for i, prot in enumerate(proteins):
        seq = prot.sequence
        seq = np.where(seq == 'U', 'C', seq) 
        L = len(seq) - k + 1
        start_range = np.arange(L)
        positions[i] = np.random.choice(start_range, p=prot.prob)   # biased sampling in intialization

    subseqs = []
    for i, prot in enumerate(proteins):
        seq = prot.sequence
        seq = np.where(seq == 'U', 'C', seq) 
        start = int(positions[i])
        std = k / 8                            # indel probability adjustable
        indel = int(np.random.normal(0, std))   # random indel (indel > 0: insert; indel < 0: delete)
        end = min(int(positions[i]) + k + indel, len(seq))
        subseq = seq[start:end]
        subseqs.append(subseq)

    # fit hmm
    model = initialize_hmm(k)
    model.fit(subseqs, transition_pseudocount=0.01, use_pseudocount=True, algorithm='viterbi')


    # gibbs sampling
    old_score = [0] * len(proteins)
    weights = np.ones(len(proteins))
    scores_100 = []     # score of every 100 iterations
    for it in range(iterations):
        # print(f'\riteration: {it}', end='')

        if it % 100 == 0:
            # score
            score = 0
            for state in model.states:
                if state.name.startswith('m'):
                    # print(f"{state.name}: {state.distribution.parameters[0]}", end=' ')
                    info = np.log2(20)
                    for r, p in state.distribution.parameters[0].items():
                        info += p * np.log2(p, where=p > 0)
                    # print('info:', info)
                    score += info   
            avg_score = score / k
            # print(f'avg score: {round(avg_score, 3)} it: {it}')
            if avg_score > 1.8:
                drop_out = True
                # print('drop out'*100)

            # calculate avg and std score (log_probability)
            if it > 1:    
                scores_100 = np.array(scores_100)
                score_avg = np.mean(scores_100[scores_100>-50])
                score_sd = np.std(scores_100[scores_100>-50])
                scores_100 = []

        i = 0 
        while i < len(proteins):
            # sample
            prot = proteins[i]
            seq = prot.sequence
            seq = np.where(seq == 'U', 'C', seq) 
            L = len(seq) - k + 1
            start_range = np.arange(L)
            start = np.random.choice(start_range, p=prot.prob)   # biased sampling in intialization
            std = k / 10                            # indel probability adjustable
            indel = int(np.random.normal(0, std))   # random indel (indel > 0: insert; indel < 0: delete)
            end = min(start + k + indel, len(seq))
            subseq = seq[start:end]

            # fit
            new_score = model.log_probability(subseq)
            old_score[i] = model.log_probability(subseqs[i])
            scores_100.append(new_score)

            # update
            if np.random.random() < np.exp((new_score - old_score[i]) / T):
                old_score[i] = new_score
                positions[i] = start
                subseqs[i] = subseq
                model = initialize_hmm(k)
                model.fit(subseqs, 
                          transition_pseudocount=0.01, 
                          use_pseudocount=True, 
                          algorithm='viterbi',
                          weights=weights)
            
            i+=1

        if it > 500:
            # drop out
            thr = score_avg + score_sd * 2      # threshold is a std above the average score (log_probability)
            for j, s in enumerate(subseqs):
                if model.log_probability(s) < thr:
                    weights[j] = 0
                else:
                    weights[j] = 1

        if it == 1000:
            # early stopping
            model_score = eval_score(model=model, k=k)
            if model_score < 1:
                break

        T *= 0.9997
    return proteins, model, weights, positions, subseqs, it


def evaluate(proteins, model, weights, positions, subseqs, k):
    evaluation = {}
    evaluation['positions'] = positions
    
    motifs = {}
    for i, prot in enumerate(proteins):
        subseq = ''.join(subseqs[i])
        motifs[prot.id] = subseq
    evaluation['motifs'] = motifs

    # score
    avg_score = eval_score(model=model, k=k)
    evaluation['avg_score'] = avg_score

    # benchmark
    benchmark = [False] * len(proteins)
    for i, prot in enumerate(proteins):
        start = positions[i]
        end = start + k
        if any(prot.motif_label[start:end]):
            benchmark[i] = True
    evaluation['benchmark'] = benchmark
    evaluation['num_benchmark'] = sum(benchmark)

    # benchmark (all specific)      benchmark of all specific motif classes 
    elms = set()
    for prot in proteins:
        for motif in prot.elm_motifs:
            elms.add(motif[0])

    benchmark_specific = {}
    for specific_motif in elms:
        d = {}
        bm_rough = 0
        bm_exact = 0
        bm_specific_tot = [False] * len(proteins)
        for i, prot in enumerate(proteins):
            for motif in prot.elm_motifs:             
                name = motif[0]
                if name == specific_motif:
                    bm_specific_tot[i] = True
                    start = motif[1]
                    if positions[i]+1 >= start - 3 and positions[i]+1 <= start + 3:
                        bm_rough += 1
                    if positions[i] == start:
                        bm_exact += 1
        bm_specific_tot = sum(bm_specific_tot)
        
        d['total'] = bm_specific_tot
        d['benchmark_rough'] = bm_rough
        d['benchmark_exact'] = bm_exact
        benchmark_specific[specific_motif] = d
    evaluation['benchmark_specific'] = benchmark_specific


    evaluation['remaining_proteins'] = int(sum(weights))
    evaluation['weights'] = weights

    return evaluation


def compute(pair, proteins, input_proteins_id, center):
    result = {}
    k, seed = pair
    seed = seed
    iterations = 3000
    k = k
    T = 10
    w_plm = 5
    w_disorder = 1
    w_rlc = 1
    w_rsa = 1

    t1 = time()
    proteins, model, weights, positions, subseqs, it = search(proteins=proteins, seed=seed,  iterations=iterations, k=k, T=T, w_plm=w_plm, w_disorder=w_disorder, w_rlc=w_rlc, w_rsa=w_rsa)
    t2 = time()
    result['seed'] = seed
    result['target_iterations'] = iterations
    result['completed_iterations'] = it
    result['k'] = k
    result['T'] = T
    result['w_plm'] = w_plm
    result['w_disorder'] = w_disorder
    result['w_rlc'] = w_rlc
    result['w_rsa'] = w_rsa
    result['input_proteins'] = input_proteins_id
    result.update(evaluate(proteins, model, weights, positions, subseqs, k=k))
    result['run_time'] = round(t2-t1, 2)

    # # To calculate computing speed
    # with open(f'{center}.txt', 'a') as f:
    #     f.write(f'{pair}\n')

    return result




