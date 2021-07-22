package service

import (
	"math/rand"
	"sort"
	"time"
)

type DNA interface {
	OrderMuta() DNA
	NodeMuta(int)

	Fitness() float64
	Len() int
	Get(int) interface{}
	Key(int) string
	Set(int, interface{})
	Copy() DNA
}

var iii int

func newDNAWithRandSort(this DNA) DNA {
	iii++
	newDna := this.Copy()
	for index := 0; index < newDna.Len(); index++ {
		newIndex := rand.Intn(newDna.Len()-index) + index
		newIndexNode := newDna.Get(newIndex)
		newDna.Set(newIndex, newDna.Get(index))
		newDna.Set(index, newIndexNode)
		//	check(newDna)
	}
	return newDna
}
func SimpleOrderMuta(dna DNA, mutaP float64) DNA {
	newDna := dna.Copy()
	for index := 0; index < newDna.Len(); index++ {
		if rand.Float64() < mutaP {
			newIndex := rand.Intn(newDna.Len()-index) + index
			newIndexNode := newDna.Get(newIndex)
			newDna.Set(newIndex, newDna.Get(index))
			newDna.Set(index, newIndexNode)
		}
	}
	return newDna
}

func (this *GA) GA(dna DNA) (best DNA) {
	endAt := time.Now().Add(this.stopTrigger.timeLimit)
	dnaList := []DNA{dna}
	var bestFitness float64 = dnaList[0].Fitness()
	best = dnaList[0]
	for len(dnaList) < this.initialPopulationNum {
		dnaList = append(dnaList, newDNAWithRandSort(dna))
	}
	this.logGen(dnaList)
	for i := 0; i < this.stopTrigger.genNum; i++ {
		//cross
		if this.stopTrigger.timeLimit > 0 && endAt.Before(time.Now()) {
			break
		}
		afterCross := this.cross(dnaList)
		this.muta(afterCross)
		better := this.keepBetterDNA(dnaList)
		bestFitness = better[0].Fitness()
		best = better[0]
		if this.stopTrigger.bestFitnesstReach != nil && bestFitness >= *this.stopTrigger.bestFitnesstReach {
			break
		}
		newGen := []DNA{}
		for len(newGen) < this.initialPopulationNum-len(better) {
			newGen = append(newGen, pickOne(afterCross))
		}
		newGen = append(newGen, better...)
		dnaList = newGen
		this.logGen(dnaList)
	}
	return best
}

type GA struct {
	initialPopulationNum int
	genealogy            genealogyStruct
	pepoleMutaP          float64
	pepoleNodeMutaP      float64
	recordGenealogy      bool
	stopTrigger          stopTrigger
	Counter              Counter
}
type genealogyStruct struct {
	dnaG            [][]DNA
	bestFitnessList []float64
}
type Counter struct {
	genNum int
	Muta   int
}
type stopTrigger struct {
	timeLimit         time.Duration
	genNum            int
	bestFitnesstReach *float64
}

func NewGA() *GA {
	return &GA{
		initialPopulationNum: 30,
		pepoleMutaP:          0.2,
		pepoleNodeMutaP:      0.2,
		stopTrigger: stopTrigger{
			timeLimit:         time.Millisecond * 1000 * 10,
			genNum:            30,
			bestFitnesstReach: nil,
		},
	}
}
func (this *GA) FitnessGenealogy() []float64 {
	return this.genealogy.bestFitnessList
}
func (this *GA) SetBestFitnessNeed(f float64) {
	this.stopTrigger.bestFitnesstReach = &f
}
func (this *GA) SetTimeLimit(l time.Duration) {
	this.stopTrigger.timeLimit = l
}
func (this *GA) SetMaxGenerNum(m int) {
	this.stopTrigger.genNum = m
}
func (this *GA) IfRecordGenealogy(r bool) {
	this.recordGenealogy = r
}
func (this *GA) logGen(dnaList []DNA) {
	if len(dnaList) < 0 {
		panic("")
	}
	this.Counter.genNum++
	if this.recordGenealogy {
		this.genealogy.dnaG = append(this.genealogy.dnaG, dnaList)
	}
	this.genealogy.bestFitnessList = append(this.genealogy.bestFitnessList, dnaList[0].Fitness())
}

type sortDNAList []DNA

func (this sortDNAList) Less(a, b int) bool {
	return this[a].Fitness() > this[b].Fitness()
}
func (this sortDNAList) Swap(a, b int) {
	this[a], this[b] = this[b], this[a]
}
func (this sortDNAList) Len() int {
	return len(this)
}
func (this *GA) keepBetterDNA(dnaList []DNA) []DNA {
	sort.Sort(sortDNAList(dnaList))
	var newGen []DNA
	for i := 0; i < 1+len(dnaList)/10 && i < len(dnaList); i++ {
		newGen = append(newGen, dnaList[i])
	}
	return newGen
}
func (this *GA) cross(dnaList []DNA) []DNA {
	sort.Sort(sortDNAList(dnaList))
	var newGen []DNA
	for {
		if len(newGen) < int(float64(len(dnaList))*1.5) {
			t := pickTwoDifferent(dnaList)
			newGen = append(newGen, crossAAndB(t[0], t[1])...)
		} else {
			break
		}
	}
	return newGen
}

func check(a DNA) {
	unionMap := make(map[string]bool)
	for i := 0; i < a.Len(); i++ {
		_, isOK := unionMap[a.Key(i)]
		if !isOK {
			unionMap[a.Key(i)] = true
		} else {
			panic("")
		}
	}
}
func crossAAndB(a, b DNA) []DNA {
	var c DNA = a.Copy()
	unionMap := make(map[string]bool)
	bitset := make([]bool, a.Len())

	for i := 0; i < a.Len(); i++ {
		if rand.Float64() > 0.5 {
			c.Set(i, a.Get(i))
			unionMap[a.Key(i)] = true
			bitset[i] = true
		}
	}
	bi := 0
	for ci := 0; ci < a.Len(); ci++ {
		if bitset[ci] {
			continue
		}
		for bi < b.Len() {
			_, isOK := unionMap[b.Key(bi)]
			if !isOK {
				unionMap[b.Key(bi)] = true
				c.Set(ci, b.Get(bi))
				bi++
				break
			} else {
				bi++
			}
		}
	}
	//	check(c)
	return []DNA{c}

}
func (this *GA) muta(dnaList []DNA) {
	sort.Sort(sortDNAList(dnaList))
	for i := 0; i < len(dnaList); i++ {
		if i < 1+len(dnaList)/5 {
			continue
		}
		if rand.Float64() < this.pepoleMutaP {
			this.Counter.Muta++
			dnaList[i] = dnaList[i].OrderMuta()
		}
	}
}

func pickOne(dnaList []DNA) DNA {
	a := dnaList[rand.Intn(len(dnaList))]
	b := dnaList[rand.Intn(len(dnaList))]
	if a.Fitness() > b.Fitness() {
		return a
	} else {
		return b
	}
}

func pickTwoDifferent(dnaList []DNA) [2]DNA {
	if len(dnaList) < 2 {
		panic("")
	}
	if len(dnaList) == 2 {
		return [2]DNA{dnaList[0], dnaList[1]}
	}
	var a, b, c int
	a = rand.Intn(len(dnaList))

	for {
		b = rand.Intn(len(dnaList))
		if b != a {
			break
		}
	}
	for {
		c = rand.Intn(len(dnaList))
		if c != a && c != b {
			break
		}
	}
	temp := []DNA{dnaList[a], dnaList[b], dnaList[c]}
	sort.Sort(sortDNAList(temp))
	return [2]DNA{temp[0], temp[1]}
}
